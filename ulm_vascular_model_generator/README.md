# DCCO ULM Vascular Tree Generator

This folder generates an X-Z planar SWC vascular tree for the ULM bubble flow
simulator using the DCCO algorithm from:

> Maso Talou et al., *Adaptive constrained constructive optimisation for
> complex vascularisation processes*, Scientific Reports 11:6180 (2021).
> https://doi.org/10.1038/s41598-021-85434-9

Default configuration file:

```text
ulm_vascular_model_generator/configs/dcco_config.yaml
```

Each run creates a timestamp-and-seed directory. For example:

```text
ulm_vascular_model_generator/vessel_swc_models/20260726_153012_123456_seed_105/
```

The directory contains the SWC file, the vessel transport NPZ file,
`radius_solver_convergence.csv`, and `radius_solver_convergence.png`. The plot
shows the fixed-point radius/viscosity error for an accepted tree update and
the number of iterations needed throughout tree growth.

The same SWC is also copied to the `vessel_swc_models` root using the configured
file name. This stable path is overwritten on each run so downstream tools can
always read the latest generated tree without knowing the timestamp directory.

The SWC uses the simulator-compatible layout:

```text
id type x y z radius parent
```

- IDs are 0-based and sequential.
- The root parent is `-1` and the root node uses SWC type `2`.
- Coordinates and radii are in micrometres.
- All nodes lie in the same X-Z plane (`y = 2000 um` by default).
- The default graph is a directed binary tree, while optional DCCO distal
  branching can represent N-ary junctions from pre-existing vessels.

## Algorithm Overview

For each new terminal, the generator:

1. Samples a candidate distal point in the perfusion domain. Uniform sampling is
   the default; Gaussian and Gaussian-mixture sampling are available for the
   paper's heterogeneous blood-flow distribution idea.
2. Rejects candidates closer than `lmin` to the current tree, with `fr`
   shrinkage after repeated failures (Eqs. 12-14).
3. Searches nearby vessels within `fn * lc`.
4. Tests candidate bifurcation points according to the vessel branching mode:
   `versatile` searches the triangle from Fig. 1, `fixed` searches the original
   vessel line, `distal` attaches at the distal endpoint, and `non_branching`
   excludes the vessel.
5. Applies Murray's law, symmetry, aspect ratio, opening-angle, domain, and
   intersection constraints (Eqs. 3-5 and 17).
6. Minimises either volume cost (Eq. 2) or optional sprouting cost (Eq. 15).

The generator has one direct flow/radius model. All terminal vessels receive an
equal share of the configured inlet flow, and segment flows are accumulated
from leaves toward the single root. Segment radii then follow the fixed-root
Murray relation directly.

The generator does not calculate viscosity, resistance, pressure, or shear.
Those properties remain the responsibility of downstream simulation stages.

## Regenerate

From the repository root:

```powershell
python .\ulm_vascular_model_generator\vessel_generation.py
```

The Python environment must provide `PyYAML` because the generator reads
parameters from YAML.

Use a custom YAML file:

```powershell
python .\ulm_vascular_model_generator\vessel_generation.py `
  --config .\ulm_vascular_model_generator\configs\dcco_config.yaml
```

To visualize one generated run, edit `RESULT_FOLDER_NAME` near the top of
`visualize_vessel_swc.py`. Set it to a folder name to pin that run, or set it
to `None` to automatically select the result folder with the latest timestamp.
Then run:

```powershell
python .\ulm_vascular_model_generator\visualize_vessel_swc.py
```

The script reads the newest SWC in that explicit run folder. Its publication
PNG/PDF/SVG, editable caption text, exported VTP mesh, and a separate PNG
matching the PyVista popup window are saved directly beside the SWC. The popup
image uses the `_pyvista_view.png` suffix and does not overwrite the flat
publication figure. By default, it also opens a resizable interactive PyVista
window: use the mouse wheel or right drag to zoom, left drag to rotate, middle
drag to pan, and `R` to reset the view. Close the window after inspection so
the script can finish. For unattended or headless rendering, add `--no-window`;
`--show-window` explicitly overrides a custom YAML file that disables it.

The default visualization configuration enables `publication_style`. The
selected `RESULT_FOLDER_NAME` is classified from
`generation_validation.json`, with the SWC coordinate span as a fallback for
older runs. A constant-Y tree is reported as a genuine two-dimensional X-Z
structure and saved as a flat orthographic figure, not as a three-dimensional
projection. Its main panel uses relative coordinates with regular
500-micrometre ticks. Every segment has the same displayed line width, while
color alone represents the physical local radius. Root, bifurcation, terminal,
and flow-direction markers define the topology. Supporting panels show a
representative bifurcation, the radius distribution, and radius against
cumulative path distance from the root.

A result with varying Y coordinates is reported as `volumetric_3d`, rendered
from an isometric three-dimensional camera, and accompanied by quantitative
X-Y, X-Z, and Y-Z projections. Both geometry modes write publication
PNG/PDF/SVG files, an editable caption text file, the PyVista screenshot, and
the VTP mesh back into the selected result folder.

The original SWC geometry and radii remain unchanged. The VTP export and
interactive PyVista cylinder view are retained for geometry inspection, but
the publication PNG is not a perspective cylinder rendering. `render.font_size`
controls annotation size, and `radius_scale` affects only the PyVista display
mesh.

All generator parameters live in the YAML file. The default
`configs/dcco_config.yaml` includes English comments explaining each variable. SWC
files are always written in a timestamp-and-seed subdirectory under
`ulm_vascular_model_generator/vessel_swc_models`; the YAML `output` value only
sets the SWC file name.

Selected parameter overrides in YAML:

```yaml
output: svv_ulm_4mm_xz_planar_dcco_tree.swc

config:
  geometry_mode: planar_2d
  seed: 15
  cube_um: 4000
  plane_y_um: 2000
  root_x_um: 0
  root_z_um: 1685
  planar_root_half_width_um: 22
  planar_min_half_width_um: 4
  planar_inlet_total_flux_um2_s: 50000
  planar_murray_gamma: 2.0
  n_terminals: 150
  delta: 0.0
  nu: 1.5
  fr: 0.9
  n_fail: 8
  delta_v: 7
  fn_neighborhood: 8.0
```

Planar and volumetric haemodynamics do not share configuration variables.
For `volumetric_3d`, use `volumetric_root_radius_um`,
`volumetric_min_radius_um`, `volumetric_inlet_total_flow_um3_s`, and
`volumetric_murray_gamma`. Planar flow is a true per-depth flux in `um^2/s`;
volumetric flow is a circular-tube volume flow in `um^3/s`.

Inlet-to-outlet main trunk generation:

```yaml
config:
  main_trunk_enabled: true
  main_trunk_endpoint_mode: random_boundary
  main_trunk_endpoint_margin_um: 180.0
  main_trunk_min_endpoint_distance_um: 2200.0
  main_trunk_axis: x
  main_trunk_position_um: null   # cube_um / 2
  main_trunk_inlet_um: 0.0
  main_trunk_outlet_um: null     # cube_um
  main_trunk_shape: meander
  main_trunk_segments: 12
  main_trunk_amplitude_um: 260.0
  main_trunk_wavelengths: 1.4
  main_trunk_noise_um: 80.0
  main_trunk_phase_rad: null
  main_trunk_branching_mode: fixed
  main_trunk_role: transport
  main_trunk_outlet_outflow: 20.0
  main_trunk_only_branching: false
```

With `main_trunk_only_branching: false`, the first side branch starts from the
main trunk and the resulting distribution tree can continue branching
recursively. Set it to `true` only for first-order side branches that attach
directly to the trunk and do not seed dense subtrees. `main_trunk_outlet_outflow`
prescribes the outlet flow carried by the trunk, matching the outlet-vessel
setup in Figure 5 of the DCCO paper.
`main_trunk_shape: meander` stores the trunk as connected SWC segments rather
than a visual-only deformation, so downstream simulation sees the curved trunk
geometry.
`main_trunk_endpoint_mode: random_boundary` samples inlet/outlet coordinates on
the domain boundary using the generator seed, while
`main_trunk_min_endpoint_distance_um` rejects pairs that are too close. Use
`fixed_axis` to recover the older deterministic x- or z-spanning trunk.

## DCCO Extensions

Staged growth is defined in the same YAML file:

```yaml
stages:
  - name: coarse_tree
    tinit: empty
    n_terminals: 50
    delta: 0.7
  - name: fine_tree
    tinit: previous
    n_terminals: 450
    delta: 0.2
```

In staged mode, each stage entry is a YAML config override, and
`n_terminals` is the number of terminals to add during that stage. Each stage
is grown on an independent clone of its `Tinit` and then unioned back into the
global tree. Every stage must define `tinit`; missing or ambiguous values raise
a configuration error.

`tinit` examples:

```yaml
stages:
  - name: first_pass
    tinit: empty
    n_terminals: 40

  - name: refine_first_pass
    tinit: "stage:first_pass"
    n_terminals: 120

  - name: grow_from_selected_outputs
    tinit:
      mode: stages
      stages: [first_pass, refine_first_pass]
      include_ancestors: true
      include_descendants: true
    n_terminals: 200

  - name: terminal_refinement
    tinit:
      mode: terminals
      terminal_ids: [12, 19, 27]
    n_terminals: 60
```

Supported `tinit` modes are `empty`, `global`, `previous`, `stage:<name>`,
explicit `vessel_ids`, explicit `terminal_ids`, and simple filters by `roles` or
`branching_modes`. Ancestors are included by default to keep the selected subset
connected; unselected ancestor connectors are temporarily non-branching inside
that stage.

Heterogeneous terminal distribution:

```yaml
config:
  terminal_distribution: gaussian
  gaussian_mean_x_um: 2000
  gaussian_mean_z_um: 2000
  gaussian_sigma_x_um: 600
  gaussian_sigma_z_um: 600
```

Avascular and carriage masks:

```yaml
config:
  avascular_regions:
    - kind: circle
      x_um: 1500
      z_um: 1500
      radius_um: 300
  carriage_regions: []
```

Multi-criteria growth:

```yaml
config:
  cost_function: sprouting
  sprout_cv: 1.0
  sprout_cp: 0.5
  sprout_cd: 1.0
```

There is no viscosity calculation or viscosity mode selector in the vascular
generator. A globally propagated pressure field is also left to the downstream
simulator.

Generated vessel role:

```yaml
config:
  default_vessel_role: perforator
```

## Branching-process GIF

Render every permanently accepted branching state using the same configuration,
sampling, constraints, flow accumulation, and staged-growth operations as the
generator:

```powershell
python .\ulm_vascular_model_generator\figure_drawing\draw_branching_process_gif.py
```

The default output is
`figure_drawing/outputs/vascular_branching_process.gif`. Use `--config` and
`--output` to select another generator configuration or destination.

## Use in the Bubble Simulator

In the MATLAB simulator, set the SWC name to the generated file stem, for
example:

```matlab
name = 'svv_ulm_4mm_xz_planar_dcco_tree';
```
