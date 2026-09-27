# Refined CFD extension surface preparation

## Scope

`cfd_surface_prepare.py` derives a CFD-specific surface from the frozen PASS preprocessing run and
the validated Ultraliser lumen. It does not rerun ROI generation, Ultraliser, the global 1D flow
solver, volume meshing, or 3D CFD. The original lumen STL and the previous direct-extrusion run are
read-only references whose SHA-256 hashes are checked before and after processing.

Only the four artificial extensions are refined. The original vessel core and the proximal cut
loop shared with that core remain exactly locked. The distal ring and cap remain locked to the
saved CFD plane. No whole-surface smoothing, remeshing, voxelization, marching cubes, implicit
reconstruction, or Ultraliser optimization is used.

## Local mesh target

For each saved boundary, original Ultraliser wall triangles are sampled within twice the source
radius on the inward side of the boundary plane. This excludes the outward rounded terminal end.
The program records local edge-length percentiles, triangle-area percentiles, and aspect-ratio
statistics in `local_original_mesh_statistics.csv`.

The local median edge length, multiplied by the YAML spacing factor, is the target axial spacing.
Each port independently receives

```text
ring_count = clamp(ceil(extension_length / target_spacing), minimum, maximum)
```

rings, including the locked proximal and distal stations. The final station is exactly the saved
extension end. Therefore a long extension can no longer be represented by only its two end rings.

## Shape transition and constrained smoothing

Ring 0 is the exact cut loop and never moves. A target loop is constructed in its two-dimensional
plane by eight low-relaxation cyclic smoothing steps. Its vertices retain cyclic order and count;
their spacing is equalized along the lightly smoothed non-circular curve. Polygon centroid and area
are restored after processing, so regularization removes mesh-scale irregularity without imposing
a perfect circle or changing vessel calibre.

Across the configured one-diameter transition, a smoothstep weight blends the actual cut loop into
the regularized target. Beyond that zone, the numerical tube retains the regularized section. The
tube remains straight, with no taper, flare, contraction, noise, voxel bumps, or artificial surface
ripple.

Six constrained Taubin iterations act only on intermediate extension rings. After every positive
and negative pass, each ring is projected back to its fixed axial station, recentered, and uniformly
scaled to its target area. The original core, Ring 0, distal ring, and cap never participate.

Every pair of neighboring rings forms quadrilateral strips. Both possible diagonals are evaluated,
and the one with the lower worst aspect ratio—then higher minimum angle—is selected. Edge collapse
is disabled. The current structured rings already satisfy the configured axial edge-size limit, so
no unsafe split of the locked interface is performed.

## Mesh quality controls

Each extension triangle is measured for its three edge lengths, area, minimum and maximum angles,
aspect ratio, and neighboring-area ratio. The reported aspect ratio is

```text
(sqrt(3) / 2) * longest_edge / altitude_to_longest_edge
```

and equals one for an equilateral triangle. A shape-bad triangle violates both the minimum-angle
and maximum-aspect criteria, while an oversized triangle independently violates the local-target
edge ratio. The aggregate bad-triangle fraction is evaluated on the refined extension body; the
locked Ring-0-to-Ring-1 interface is excluded from that aggregate and subjected to its own explicit
edge-ratio, P95 aspect, and minimum-angle report. Neighbor-area P95 is checked separately. This
separation prevents immutable microscopic Ring 0 edges from being hidden while avoiding double
counting the independently reviewed interface.

The run also requires improvement over the previous direct extrusion: every port must have more
triangles, lower P95 aspect ratio, and a higher minimum angle. Existing topology, collision, dense
core closest-point, exact original-vertex motion, cap geometry, source hashes, and meter-scale
checks remain mandatory.

## Boundary conditions

The inlet remains the saved parabolic volumetric-flow boundary. Each outlet pressure is recalculated
from the final refined cap area—not copied from the old run—using only the Poiseuille loss of the
artificial extension:

```text
R_ext = 8 mu L_ext / (pi r_eq^4)
delta_P_ext = R_ext Q_1D
P_solver = P_original - delta_P_ext
```

Original P/Q values remain unchanged and traceable. This is numerical artificial-extension
correction, not a physiological outlet model, resistance boundary condition, or Windkessel model.
Negative gauge pressure at the true terminal remains valid.

## Command and outputs

Run from the project directory after tests and Ruff have passed:

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_surface_prepare.py
```

The strict default configuration is `configs/cfd_surface_prepare.yaml`. A successful run creates:

```text
outputs/cfd_surface_prepare/cfd_surface_refined_anchor003274_<timestamp>/
  input/
    cfd_surface_prepare.yaml
    source_manifest.json
    original_surface_reference.json
    cfd_preprocess_reference.json
    previous_direct_extrusion_reference.json
  geometry/
    cfd_surface_refined_um.stl
    cfd_surface_refined_um.vtp
    cfd_surface_refined_m.stl
    cfd_wall_refined_um.stl
  boundaries/
    boundary_00_*.stl ... boundary_03_*.stl
    boundary_manifest.csv
  bc/
    boundary_conditions_original.json
    boundary_conditions_refined.json
    extension_pressure_correction_refined.csv
  qc/
    local_original_mesh_statistics.csv
    extension_mesh_quality_qc.json
    extension_mesh_before_after.csv
    original_locked_vertex_motion_qc.json
    boundary_geometry_qc.json
    core_surface_preservation_qc.json
    surface_qc.json
    extension_collision_qc.json
    original_surface_integrity.json
    meter_scale_qc.json
    run_summary.json
  figures/
    previous_vs_refined_surface.png
    extension_interface_before_after.png
    extension_wireframe_before_after.png
  report/
    cfd_surface_prepare_report.md
```

The VTP retains wall/inlet/outlet, boundary index, boundary origin, extension index, and axial-band
cell tags. The four distal cap STLs remain separate boundary patches. A successful automatic result
is deliberately named `CFD_EXTENSION_MESH_REFINED_PENDING_MANUAL_REVIEW`; it is not a declaration
of volume-mesh or CFD-solver readiness. Processing stops at
`MANUALLY REVIEW REFINED CFD STL SURFACE`.
