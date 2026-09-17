# Microbubble Flow Reproduction Specification

This file documents the current `ulm_microbubble_traj_gen` implementation.
It is a reproduction checklist for the active field-based pipeline: inputs,
configuration, numerical model, diagnostics, outputs, and commands required to
repeat a run.

Implementation state documented here: **2026-07-19**.

This revision was audited against the executable entry point, its current
default YAML, the complete `utils` call chain, output writers,
visualization loaders, and the tests present in the working tree. Values called
"shipped YAML" below are deliberately distinguished from loader fallbacks used
when a field is omitted.

## 0. Current Scope

| Item                    | Current implementation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Main task               | Generate finite-size hydrodynamic-mobility microbubble trajectories from a formally empty lumen using deterministic concentration-driven continuous perfusion.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Active geometry         | 2D X-Z lumen mask derived from the DCCO vessel graph. The saved particle Y coordinate is fixed.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| Active flow model       | Hybrid PhiFlow/SciPy 2D finite-volume lumen flow. PhiFlow performs explicit viscous relaxation; SciPy CG performs face-flux pressure projection.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Active cardiac model    | The accepted steady CFD field is the cycle-mean reference. A positive ECG-shaped surrogate is mean-normalized, evaluated continuously at every internal particle stage, and retarded by root-path distance divided by pulse-propagation speed. This is a quasi-steady kinematic modulation, not compliant transient CFD.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Active solver mode      | `phiflow_viscous_fv_projection_2d`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| Active particle model   | Position-, radius-, viscosity-, and exact-wall-gap-dependent mobility with near-wall shear forcing, deterministic pair collisions, synchronous event-split integration, and the Revised-v6 equal-flux schedule. Revised v15 defines one lumen domain\(\Omega\), one finite-segment solid-wall set \(\Gamma_w\) that excludes the authoritative CFD open faces, and directed terminal outlet sections \(\psi_k\). Finite-size feasibility uses \(g_R=d_w-R\), while hydrodynamic regularization is confined to mobility coefficients. Every internal move predicts penetration before movement and applies the closed-form single-wall mobility reaction; a non-feasible complete chord triggers genuine physical-time bisection. Outlet crossing has first priority, continuing particles must remain in the inlet-connected finite-radius domain \(\Omega_R^{\mathrm{in}}\), and a true simultaneous two-wall contact is an explicit model-limit error. The optional molecular-target and deterministic mean-field bond extensions remain active as configured. Formal frame 0 is empty. |
| Removed particle mode   | `passive_tracer` and the former fixed-population transport are no longer present. Configuration and public dispatch expose only concentration-driven continuous perfusion.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Active particle backend | NumPy lifecycle/event orchestration with fused Numba kernels for hydrodynamic mobility, applied loads, collision pairs, optional molecular state, exact finite-face distance, and swept-disc wall audits. The recorded kernel family is`numba_batched_component_kernels_v18`; deterministic outlet ordering and physical-time refinement remain orchestration-level. The `python` backend executes the same public model as a reference path. Taichi is not selected for the current synchronous CPU batches.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| Required vessel input   | `.vessels.npz` exported by `ulm_vascular_model_generator/vessel_generation.py`; `.swc` is loaded as a companion path.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Hard stop conditions    | Missing dependencies, invalid configuration, disconnected/holed lumen, nonconverged flow, invalid canonical boundary data, invalid inlet section or zero MB flux, deterministic-ID/record limits, invalid mobility/collision numerics, a continuing centre outside\(\Omega\) or \(\Omega_R^{\mathrm{in}}\), unresolved position complementarity, true simultaneous distinct-wall contact, or incomplete formal streamline coverage.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| Not active now          | Legacy diffusion-only pseudo-flow, thick inlet/outlet velocity bands, post-projection flow cleanup, the Revised-v13 velocity-projection/closest-point-repair path, persistent wall-contact active sets, arbitrary nearest-wall selection at simultaneous two-wall contact, displacement clipping, zero-progress acceptance, discrete stochastic bond identities, target depletion/internalization, acoustic/buoyancy states or forces, long-range many-body hydrodynamic coupling, full 3D flow, and MAT output.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

The active command-line entry point is:

```text
ulm_microbubble_traj_gen/generate_microbubble_trajectories.py
```

## 1. End-to-End Workflow

```text
parse --config, --quick-test, --skip-render, optional --reuse-field-from,
and optional --prepare-target-candidates
    -> load and validate YAML, resolve paths, and construct the timestamped output path
    -> validate PyVista/VTK/trame visualization dependencies
    -> create the output directory
    -> optionally write the source-derived run_config.yaml
    -> exported DCCO vessel model
    -> load .swc and .vessels.npz
    -> build an X-Z grid domain
    -> generate a continuous lumen polygon from DCCO centerlines and radii
    -> compute exact Shapely cell/lumen intersection fractions and a centre-valid majority mask
    -> assign vessel attributes to lumen cells without changing the mask
    -> require one connected lumen component and holes=0
    -> build root-inlet and terminal-outlet open-edge face fluxes
    -> iterate PhiFlow viscous relaxation + SciPy finite-volume projection
    -> accept only a physically converged velocity field
       -> on nonconvergence: optionally save the failed field, always write
          failure diagnostics and domain metadata, then abort
    -> write diagnostics for the accepted field before particle transport
    -> either:
       -> candidate-preparation mode: build flow-oriented vessel-bed units,
          candidate masks and metrics, save the compact catalog and final VTI,
          then stop before particle transport
       -> formal simulation mode: precompute gradients, viscosity extension,
          solid-wall distance, and inward normal
    -> build one canonical particle boundary geometry from the authoritative CFD open faces:
       lumen Omega, solid faces Gamma_w, directed outlet sections psi_k, and Omega_R^in
    -> optionally intersect a coordinate-aware Boolean mask with Gamma_w
    -> build one periodic cardiac-flow cycle and a root-path propagation-delay field
    -> construct concentration-calibrated inlet flux and integrated pulsatile injection events
    -> define formal frame 0 as an empty lumen
    -> synchronously evaluate background hydrodynamics + pair collisions + optional mean-field molecular bonds
    -> predict and solve each wall reaction with the Revised-v15 closed-form mobility constraint
       using true gap g_R; use the regularized hydrodynamic gap only inside mobility coefficients
    -> compare directed outlet and wall events, give the outlet an exact-tie priority,
       and terminate the permanent ID exactly on the outlet plane
    -> bisect failed complete chords into two real half-time intervals, recomputing all stage physics
    -> require every continuing centre to remain in Omega and Omega_R^in
    -> optionally report no-bond target exposure and predeclared Da_on scenarios whose
       reference time was fixed in configuration before transport
    -> conditionally save field, trajectory, and target NPZ files
    -> always write successful-run domain metadata
    -> render the final WSS scene
    -> render initial/final CFD scenes and validate formal root-to-outlet streamline coverage
    -> return result paths and print completion/viewer commands
```

| Step                              | Input                                                                                                                  | Output                                                                                                                                                                                                                                                        | Main code                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| --------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1. Export vascular model          | DCCO vascular-generator config                                                                                         | `<model>.swc`, `<model>.vessels.npz`                                                                                                                                                                                                                      | `ulm_vascular_model_generator/vessel_generation.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| 2. Load vessels                   | `.swc`, `.vessels.npz`                                                                                             | Vessel objects with geometry, radius, flow, viscosity, pressures, parent/children                                                                                                                                                                             | `utils/vascular_io.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| 3. Build grid domain              | Vessel coordinates,`domain` config                                                                                   | X-Z grid origin, spacing, shape, fixed Y                                                                                                                                                                                                                      | `utils/grid_domain.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| 4. Rasterize lumen                | Vessel centerlines, radii, topology                                                                                    | `lumen_mask`, `wall_mask`, `junction_core_mask`, vessel IDs, per-cell physical fields                                                                                                                                                                   | `utils/vessel_rasterizer.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| 5. Check mask quality             | Rasterized lumen mask                                                                                                  | Error unless one connected component and`holes=0`                                                                                                                                                                                                           | `utils/connectivity.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| 6. Build open boundaries          | Mask boundary faces, root and terminal vessels                                                                         | Per-face inlet/outlet flux labels, normals, lengths, target fluxes                                                                                                                                                                                            | `utils/flow_boundaries.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| 7. Solve flow                     | Raster, boundaries, field config                                                                                       | Velocity, speed, pressure, face fluxes, divergence, wall penetration, WSS                                                                                                                                                                                     | `utils/phiflow_solver.py`, `utils/face_flux_projection.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| 8. Diagnose flow                  | Accepted or failed flow                                                                                                | YAML summary, divergence heatmap, outlet flux CSV, narrow segment CSV                                                                                                                                                                                         | `utils/flow_diagnostics.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| 9. Prepare particle hydrodynamics | Accepted grid velocity/viscosity plus the pre-raster continuous vessel geometry                                               | Grid velocity gradient and extended Pa-s viscosity for interpolation; one required continuous geometry for wall distance, normal, contact, inlet/outlet, and finite-radius accessibility                                                                 | `utils/particle_hydrodynamic_fields.py`, `utils/continuous_vessel_geometry.py`, `utils/particle_boundary_support.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| 9A. Candidate target preparation  | Accepted CFD, vessel topology, raster ownership, solid wall, finite-size inlet flux, formal observation duration       | Flow-oriented vessel-bed units, physical wall-area quadrature, local/subtree candidates, accessibility metrics, optional compact influence region and correlated wall patches, compact NPZ/JSON catalog, and selector VTI                                     | `utils/vessel_bed_topology.py`, `utils/molecular_target_wall_measure.py`, `utils/molecular_target_candidates.py`, `utils/molecular_target_auto_selection.py`, `utils/molecular_target_spatial_heterogeneity.py`, `utils/molecular_target_candidate_io.py`, `utils/target_candidate_runner.py`                                                                                                                                                                                                                                                                          |
| 9B. Optional molecular target     | Solid-wall geometry plus a coordinate-aware Boolean NPZ mask                                                           | Eligible target-positive wall sites, density field, physical coordinates, and inward normals                                                                                                                                                                  | `utils/molecular_target_field.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| 10. Advance bubbles               | Accepted grid flow field, continuous boundary geometry, cardiac config, particle/dynamics config, optional target/binding config | Permanent-ID trajectories, directed outlet lifecycle, continuous-geometry predictive reaction, rotation, collision/contact/kinematic diagnostics, and optional bond states and loads                                                                      | `cardiac_waveform.py`, `utils/cardiac_pulsatility.py`, `utils/particle_inlet_flux.py`, `utils/particle_perfusion_schedule.py`, `utils/particle_perfusion_transport.py`, `utils/particle_mobility_transport.py`, `utils/particle_mobility.py`, `utils/particle_predictive_contact.py`, `utils/continuous_vessel_geometry.py`, `utils/particle_boundary_numba.py`, `utils/particle_constrained_step.py`, `utils/particle_kinematic_diagnostics.py`, `utils/particle_collisions.py`, `utils/molecular_binding.py` |
| 10A. Optional contact pilot       | Empty-start no-bond trajectories, target field, predeclared`Da_on` axes, and fixed `da_on_reference_time_s`        | Exposure report plus scenario definitions independent of observed contact time; no fitted biological rate                                                                                                                                                     | `utils/molecular_contact_pilot.py`, `utils/molecular_binding_scenarios.py`, `utils/molecular_pilot_runner.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| 11. Save results                  | Domain, raster, flow, trajectories                                                                                     | Compressed NPZ files and metadata YAML                                                                                                                                                                                                                        | `utils/field_io.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| 12. Render native CFD             | Domain, raster, initial/final flow                                                                                     | Final WSS HTML/PNG plus initial/final CFD HTML/PNG, VTI/VTP, and continuity CSV                                                                                                                                                                               | `vis_utils/pyvista_wall_shear.py`, `vis_utils/pyvista_flow.py`, `vis_utils/vtk_streamlines.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| 13. Optional trajectory rendering | Saved field/trajectory NPZ files                                                                                       | Snapshot, optional GIF/MP4, field maps, holes map, and narrow-cells map                                                                                                                                                                                       | `visualize_microbubble_flow.py`, `vis_utils/cli.py`, `vis_utils/plotting.py`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |

The microbubble generator does not run the vascular generator internally. Run
the vascular generator first when the required `.vessels.npz` file is missing
or stale.

At the command-line boundary, `ImportError` and `ValueError` are reported as
`Generation stopped: ...` and exit status 1 without a traceback.
`FlowConvergenceError` derives from `ValueError`, so it follows this path after
its diagnostic files are written. Other uncaught exception classes retain a
normal traceback. Completion text and viewer commands appear only after the
entire runner, including mandatory rendering, succeeds.

## 2. Inputs

| Input                     | Meaning                                                                                                                                                                                                                              | Default source                                                                                                                             |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `.swc`                  | Human-readable centerline companion file. It is useful for inspection but is not the complete physics input.                                                                                                                         | `ulm_vascular_model_generator/vessel_swc_models/<model>.swc`                                                                             |
| `.vessels.npz`          | Machine-readable vessel export used by the physics pipeline. This is the required model file.                                                                                                                                        | `ulm_vascular_model_generator/vessel_swc_models/<model>.vessels.npz`                                                                     |
| Physics YAML              | Selects model name, grid resolution, mask controls, flow solver tolerances, particles, outputs, and quick-test overrides.                                                                                                            | `ulm_microbubble_traj_gen/configs/physics_flow_config.yaml`                                                                              |
| Optional target mask NPZ  | Coordinate-aware Boolean target ROI used when`molecular_target.region_mode=mask_npz`; it must carry strictly increasing physical X/Z axes in micrometers.                                                                          | Prefer the Revised-v8/v10 selector or automatic-preparation output; externally justified masks remain valid when they use the same schema. |
| Optional scenario overlay | Literature-provenance or dimensionless-sensitivity fragments. These are not stand-alone production configurations and must be merged into a full YAML only after a real target ROI and dimensionally compatible values are supplied. | `ulm_microbubble_traj_gen/configs/molecular_binding_scenarios/`                                                                          |
| Optional reusable field   | A previously accepted `velocity_and_wall_shear_field.npz`, supplied directly or through its result directory. Reuse is accepted only after exact configuration, domain, raster, continuous-geometry hash, convergence, and array-contract checks. | `ulm_microbubble_traj_gen/results/<accepted-run>/`                                                                                  |

The current pipeline uses exported vessel fields directly. It does not infer
flow rate, viscosity, radius, or topology from the SWC file.

The current generator validates the CFD rendering stack for its default,
fully rendered run. Passing `--skip-render` bypasses that visualization-only
validation and rendering stage while retaining every numerical output. NumPy, SciPy, PhiFlow, Shapely,
PyYAML, and Matplotlib are required by the numerical and diagnostic paths;
Matplotlib is used even for failed-flow diagnostics. Numba accelerates supported
kernels but has a Python/NumPy fallback, and `tqdm` has a no-op fallback. The
startup visualization check requires PyVista `>=0.48`, the APIs used by this
renderer, VTK with the OpenGL2 LIC module, trame, trame-vtk, and
trame-vuetify. The reproducible visualization requirements file additionally
pins PyVista below `0.49`; that upper bound is an environment pin, not a
separate runtime version check:

```text
ulm_microbubble_traj_gen/requirements-cfd-visualization.txt
```

Install it into the project environment with:

```powershell
D:\anaconda3\envs\pmp\python.exe -m pip install -r ulm_microbubble_traj_gen\requirements-cfd-visualization.txt
```

That file is only the CFD-visualization supplement, not a complete bootstrap
file for an empty environment. GIF output also needs Pillow, while MP4 output
needs a system-accessible FFmpeg writer. The implementation was verified in
the following environment snapshot; these are observed versions rather than a
claim that every package must be hard-pinned to exactly this version:

| Component                         | Verified version                                   |
| --------------------------------- | -------------------------------------------------- |
| Python executable                 | `D:\anaconda3\envs\pmp\python.exe` (`3.10.18`) |
| NumPy / SciPy                     | `2.2.6` / `1.15.3`                             |
| PhiFlow / Shapely / PyYAML        | `3.4.0` / `2.1.2` / `6.0.2`                  |
| Matplotlib / Numba / tqdm         | `3.10.5` / `0.61.2` / `4.67.1`               |
| PyVista / VTK                     | `0.48.4` / `9.6.2`                             |
| trame / trame-vtk / trame-vuetify | `3.13.2` / `2.11.8` / `3.2.2`                |
| Pillow / packaging                | `11.3.0` / `25.0`                              |

The production `auto` particle backend selects deterministic Numba CPU kernels.
Taichi 1.7.4/CUDA was evaluated on the target machine, but it is not selected
for these small (typically tens of particles), event-driven, float64 contact
batches: per-substep GPU launch and transfer costs exceed the useful work and
would weaken deterministic event ordering. A future Taichi backend would need
device-resident state and multi-substep time blocks rather than one launch per
internal step.

Passing `--reuse-field-from <result-directory-or-npz>` skips the CFD solve only.
The loader requires the sibling source `run_config.yaml`, rebuilds the current
domain/raster/continuous geometry, and rejects any mismatch before particle
transport starts. Diagnostics, hydrodynamic preparation, transport, saving,
and optional rendering still run normally. Omitting the option preserves the
ordinary solve path.

## 3. Outputs

With the shipped YAML, runtime outputs are written under:

```text
ulm_microbubble_traj_gen/results/<timestamp>/
```

A custom absolute `output.results_dir` is used directly. A relative value is
resolved against the `ulm_microbubble_traj_gen` package directory.

| Output file or folder                                | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| ---------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `velocity_and_wall_shear_field.npz`                | Accepted grid fields plus the`v15_exact_solid_face_particle_boundary` schema: domain, masks, vessel attributes, velocity, speed, WSS, pressure, divergence, face fluxes, labels, authoritative solid/open face geometry, target/actual fluxes, and solver metadata.                                                                                                                                                                    |
| `failed_velocity_and_wall_shear_field.npz`         | Written only when the velocity field fails physical convergence but`output.save_npz` is enabled. Particle advection is not run in this case.                                                                                                                                                                                                                                                                                           |
| `microbubble_field_trajectories.npz`               | Append-only permanent-ID observations and lifecycle registry using schema v9 without molecular binding or v10 with binding; stores last accepted internal v15 constrained velocity, realized centre velocity, exact true gap, v15 contact/lifecycle diagnostics, per-frame population, metadata, and optional mobility/rotation/collision/bond arrays.                                                                                   |
| `molecular_target_candidates.npz`                  | Candidate-preparation-only schema-v3 catalog: physical axes, masks, unit ownership/topology, ragged candidate memberships, physical wall-area weights, wall-to-segment mapping, accessibility mask, and candidate metrics. It deliberately does not store one dense mask per candidate.                                                                                                                                                  |
| `molecular_target_candidates.json`                 | Human-readable companion index containing the unit hierarchy, candidate labels/memberships, and descriptive flow metrics.                                                                                                                                                                                                                                                                                                                |
| `selected_molecular_target_mask.npz`               | Manual-selector or configured automatic-selection output containing`x_um`, `z_um`, and the final Boolean `target_mask`. Manual output retains simplified candidate IDs; automatic output stores the influence masks, anchor, physical radius, requested/achieved wall fractions, correlation settings, random coefficients, field samples, threshold, and patch count. It becomes the formal `mask_npz` input only after review. |
| `automatic_molecular_target_selection.json`        | Revised-v10 automatic-mode audit report containing coarse-anchor eligibility, physical influence-region metrics, target-positive patch area, correlation realization, analytical accessibility context, mapped/unmapped wall area, and the explicit exclusion of speed, WSS, and residence time from target-positive probability.                                                                                                        |
| `molecular_target_field.npz`                       | Written during a formal run when a molecular target is enabled and`output.save_npz=true`; records the grid overlay, density, normalized source mask, and canonical \(\Gamma_w\) face centres, normals, axes, and lengths used for target reaction geometry.                                                                                                                                                                            |
| `molecular_contact_pilot_*_capture_ratio_*.yaml`   | Optional no-bond exposure report and predeclared dimensionless-`Da_on` scenario table, one per capture/rest-length ratio. Every table uses the configured pre-run `da_on_reference_time_s`, not observed contact time.                                                                                                                                                                                                               |
| `domain_metadata.yaml`                             | Human-readable run summary: grid size, spacing, input files, particle model/schema/velocity semantics, population and wall metrics, solver mode, solver metadata, and failed-field path when applicable.                                                                                                                                                                                                                                 |
| `run_config.yaml`                                  | A deep copy of the source YAML mapping when`output.save_run_config` is true, augmented with `_resolved_source_config_path` and a resolved strict field-reuse contract; a configured molecular mask path is replaced by its resolved absolute path.                                                                                                                                    |
| `flow_diagnostics/flow_diagnostics.yaml`           | Scalar diagnostics for divergence, mask quality, fluxes, solver metadata, and generated diagnostic files.                                                                                                                                                                                                                                                                                                                                |
| `flow_diagnostics/divergence_heatmap_overlay.png`  | Divergence heatmap over segment ID, lumen, inlet, outlet, and wall masks.                                                                                                                                                                                                                                                                                                                                                                |
| `flow_diagnostics/outlet_flux_errors.csv`          | Per-outlet target flux, actual final flux, and relative error.                                                                                                                                                                                                                                                                                                                                                                           |
| `flow_diagnostics/narrow_segments.csv`             | Per-segment narrow-cell statistics and diameter diagnostics.                                                                                                                                                                                                                                                                                                                                                                             |
| `initial_flow_field_cfd.html/.png`                 | Offline interactive scene and static preview for the initial velocity field and zero pressure reference.                                                                                                                                                                                                                                                                                                                                 |
| `final_flow_field_cfd.html/.png`                   | Offline interactive scene and static preview for the accepted velocity and projection-pressure fields.                                                                                                                                                                                                                                                                                                                                   |
| `initial/final_flow_field.vti`                     | Native VTK image data for quantitative field inspection. Initial WSS is unavailable and stored as NaN; final VTI contains the accepted-velocity in-plane WSS proxy.                                                                                                                                                                                                                                                                      |
| `initial/final_streamlines.vtp`                    | Complete root-to-outlet streamlines that passed continuity and topology validation.                                                                                                                                                                                                                                                                                                                                                      |
| `initial/final_streamlines_diagnostic.vtp`         | Rejected forward candidates retained only for diagnosis.                                                                                                                                                                                                                                                                                                                                                                                 |
| `initial/final_root_to_outlet_continuity.csv`      | Human-readable streamline continuity, termination, and vessel-coverage audit.                                                                                                                                                                                                                                                                                                                                                            |
| `final_wall_shear_stress_cfd.html/.png`            | Final closed-wall-band visualization of the 2D in-plane wall-shear-stress proxy in Pa.                                                                                                                                                                                                                                                                                                                                                   |
| `lumen_holes.png`                                  | High-contrast static map of holes in the saved lumen mask.                                                                                                                                                                                                                                                                                                                                                                               |
| `narrow_lumen_cells.png`                           | High-contrast static map of current narrow-cell regions.                                                                                                                                                                                                                                                                                                                                                                                 |
| `initial_flow_field.png`, `final_flow_field.png` | Matplotlib field maps produced by`visualize_microbubble_flow.py`. These are separate from the PyVista CFD previews.                                                                                                                                                                                                                                                                                                                    |
| `microbubble_flow_snapshot.png`                    | Optional static trajectory visualization.                                                                                                                                                                                                                                                                                                                                                                                                |
| `microbubble_flow.gif` or `microbubble_flow.mp4` | Optional animation output.                                                                                                                                                                                                                                                                                                                                                                                                               |

CFD visualization dependencies are validated before the result directory and
flow solve are started. For an accepted flow, diagnostics are written first,
particle transport runs, and an enabled molecular contact pilot runs **before**
the core accepted-field/trajectory NPZ files are saved. Successful-run domain
metadata is then written before final WSS and initial/final CFD rendering. A
late streamline-continuity or rendering failure can therefore leave valid core
NPZ files in the timestamp folder while the generator command still exits as a
failure. A target, particle, or pilot failure happens earlier and can leave
accepted-flow diagnostics without those core NPZ files.

When `output.save_npz=false`, the returned field and trajectory paths are only
expected locations and the CLI still prints them; those files are not created.
The returned molecular-target path is instead `None` and is not printed.
Metadata, contact-pilot YAML, and rendering remain enabled independently of
`save_npz`.

| Failure stage                                             | Files that may already exist                                                                                                                                             |
| --------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| YAML read, conversion, or configuration validation        | No timestamp directory is created.                                                                                                                                       |
| PyVista/VTK/trame dependency validation                   | No timestamp directory is created.                                                                                                                                       |
| Vessel loading, rasterization, or connectivity validation | Timestamp directory and, when enabled,`run_config.yaml`.                                                                                                               |
| Flow physical nonconvergence                              | Prior files, flow diagnostics, and`domain_metadata.yaml`; the failed-field NPZ exists only when `save_npz=true`. Metadata still records its expected path otherwise. |
| Target build, particle transport, or contact pilot        | Prior files and accepted-flow diagnostics; successful core NPZ and successful-run metadata have not yet been written.                                                    |
| Core save or metadata write                               | A partial set of core files may exist.                                                                                                                                   |
| WSS, streamline-continuity, or CFD rendering              | Core NPZ files and successful-run metadata already exist; partial native/rendering artifacts may also exist, but the command still fails.                                |
| Complete success                                          | The entry point prints completion, output paths, and viewer commands only after`run_generation()` returns.                                                             |

## 4. Configuration

Default configuration:

```text
ulm_microbubble_traj_gen/configs/physics_flow_config.yaml
```

| YAML section                   | Field                                                                                                                                            | Current meaning                                                                                                                                                                                                                                                                                                                     |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `input`                      | `model_name`                                                                                                                                   | Base name used to find`<model>.swc` and `<model>.vessels.npz`.                                                                                                                                                                                                                                                                  |
| `input`                      | `swc_dir`                                                                                                                                      | Directory that contains the exported vascular model files.                                                                                                                                                                                                                                                                          |
| `output`                     | `results_dir`                                                                                                                                  | Parent directory for timestamped result folders.                                                                                                                                                                                                                                                                                    |
| `output`                     | `timestamp_format`                                                                                                                             | Timestamp naming format for each result folder.                                                                                                                                                                                                                                                                                     |
| `output`                     | `save_run_config`, `save_npz`                                                                                                                | Enable copied YAML and compressed NPZ output.                                                                                                                                                                                                                                                                                       |
| `domain`                     | `grid_spacing_um`                                                                                                                              | Main grid resolution. The current default is`1.5 um`.                                                                                                                                                                                                                                                                             |
| `domain`                     | `padding_um`                                                                                                                                   | Extra margin around the vascular model before rasterization.                                                                                                                                                                                                                                                                        |
| `domain`                     | `min_lumen_radius_cells`                                                                                                                       | Minimum geometric rasterization radius in grid cells. It is a connectivity fallback, not a replacement for DCCO physical radius.                                                                                                                                                                                                    |
| `domain`                     | `min_resolved_diameter_cells`                                                                                                                  | Resolution warning threshold. The current default is`10.0` cells.                                                                                                                                                                                                                                                                 |
| `domain`                     | `max_grid_cells`                                                                                                                               | Safety limit that stops overly large grids.                                                                                                                                                                                                                                                                                         |
| `field`                      | `solver_mode`                                                                                                                                  | Canonical value:`phiflow_viscous_fv_projection_2d`. The implementation rejects the misleading `phiflow_incompressible_2d` name.                                                                                                                                                                                                 |
| `field`                      | `effective_thickness_um`                                                                                                                       | Converts exported 3D flow`Q` to 2D flux target `Q / effective_thickness_um`.                                                                                                                                                                                                                                                    |
| `field`                      | `boundary_depth_cells`                                                                                                                         | Search depth near root and terminal endpoints for open boundary faces.                                                                                                                                                                                                                                                              |
| `field`                      | `max_iterations`, `min_iterations`                                                                                                           | Iteration bounds for viscous relaxation plus projection.                                                                                                                                                                                                                                                                            |
| `field`                      | `divergence_tolerance`                                                                                                                         | Limit for normalized final and interior finite-volume divergence.                                                                                                                                                                                                                                                                   |
| `field`                      | `flux_tolerance`                                                                                                                               | Limit for total, target, inlet-label, and outlet-label final flux errors.                                                                                                                                                                                                                                                           |
| `field`                      | `max_wall_penetration_um_s`                                                                                                                    | Maximum allowed normal velocity on solid wall faces.                                                                                                                                                                                                                                                                                |
| `field`                      | `momentum_residual_tolerance`                                                                                                                  | Limit for the normalized discrete viscous / pressure-projection momentum-balance residual on the interior lumen.                                                                                                                                                                                                                   |
| `field`                      | `kinematic_viscosity_um2_s`                                                                                                                    | Diffusivity used by PhiFlow explicit viscous relaxation.                                                                                                                                                                                                                                                                            |
| `field`                      | `viscous_relaxation_dt_s`                                                                                                                      | Time step for each viscous relaxation pass.                                                                                                                                                                                                                                                                                         |
| `field`                      | `pressure_solver_tolerance`, `pressure_solver_max_iterations`                                                                                | SciPy CG pressure projection controls.                                                                                                                                                                                                                                                                                              |
| `particles`                  | `inlet_number_concentration_mb_per_ml`                                                                                                         | Physical inlet number concentration. The current checked-in YAML contains`2e7 MB/mL`; the commented reference conditions `5e5` and `1e6 MB/mL` are not active values. Internally the configured value is multiplied by `1e6` to obtain `bubble/m^3`.                                                                      |
| `particles`                  | `n_steps`, `dt_s`                                                                                                                            | Number and duration of formal saved-frame intervals. Defaults`1000` and `0.001 s`, producing `1001` formal frames with an empty frame 0 followed by concentration-driven filling.                                                                                                                                             |
| `particles`                  | `bubble_diameter_min_um`, `bubble_diameter_max_um`                                                                                           | Required bounds of the uniform**prior** diameter distribution, currently `1.5..2.5 um`. The public schedule samples the marginal finite-size number flux, so accessibility and inlet speed can make scheduled diameters non-uniform within those bounds. Each permanent ID keeps its assigned diameter.                     |
| `particles`                  | `contact_geometry_tolerance_um`                                                                                                                | Revised-v15 upper bound for a small exact-gradient endpoint residual projection. The checked-in production YAML uses`0.0015 um`; the loader fallback is `0.001 um`. It never authorizes an accepted negative true gap, creates a physical stand-off, or replaces the independent hydrodynamic regularization gap.               |
| `particles`                  | `contact_max_time_refinements`                                                                                                                 | Maximum recursive bisection depth for a physical-time interval whose curved-contact trial cannot yet be accepted. Each retry uses two genuine half-time intervals whose durations sum to the parent interval; exhaustion is an explicit numerical error, never an accepted zero-progress state.                                     |
| `particles`                  | `wall_contact_threshold_um`                                                                                                                    | Output-only threshold for labelling a stored observation as close to the wall. Default`0.05 um`; it neither activates the unilateral reaction nor changes geometric feasibility.                                                                                                                                                  |
| `particles`                  | `max_unique_bubbles`                                                                                                                           | Optional permanent-ID safety cap.`0` derives capacity from the deterministic schedule over the formal duration. Terminated slots are never reused.                                                                                                                                                                                |
| `particles`                  | `max_particle_frame_records`                                                                                                                   | Runtime guard on the sum of the true variable active population across formal frames. Current default`5,000,000`.                                                                                                                                                                                                                 |
| `particles`                  | `acceleration_backend`                                                                                                                         | `auto` uses installed Numba CPU kernels and otherwise falls back to Python/NumPy. The explicit values are `numba_cpu` and `python`.                                                                                                                                                                                           |
| `particle_dynamics`          | `time_integrator`                                                                                                                              | `euler` or strict `heun`. Heun recomputes fields, mobility, collisions, and the unilateral reaction at its feasible predictor state. Default `euler`.                                                                                                                                                                         |
| `particle_dynamics`          | `integration_substeps`                                                                                                                         | Fixed nominal mobility subintervals between two saved frames. The shipped YAML uses`80`, so `dt_s=0.001 s` gives `internal_dt_s=1.25e-5 s`. A field omitted from another YAML falls back to `1`; that fallback is not the shipped setting.                                                                                  |
| `particle_dynamics`          | `near_wall_enabled`                                                                                                                            | Enables wall/free-space mobility blending and wall-shear background forcing. When false, background translation is exactly the sampled flow velocity.                                                                                                                                                                               |
| `particle_dynamics`          | `collisions_enabled`                                                                                                                           | Enables deterministic mobility-scaled pairwise soft repulsion.                                                                                                                                                                                                                                                                      |
| `particle_dynamics`          | `xi_min`, `xi_near`, `xi_far`                                                                                                              | Engineering bounds for logarithm regularization, near-wall asymptote freezing, and smooth transition to free space. Defaults`1e-3`, `0.1`, `1.0`.                                                                                                                                                                             |
| `particle_dynamics`          | `collision_layer_um`                                                                                                                           | Extra short-range collision layer added to`R_i + R_j`. Default `0.05 um`; this is an engineering parameter, not a value reported by the reference paper.                                                                                                                                                                        |
| `particle_dynamics`          | `collision_relaxation_time_s`                                                                                                                  | Overlap-relaxation time used to derive pair stiffness from relative mobility. Default`0.01 s`; keep fixed during time-step convergence studies.                                                                                                                                                                                   |
| `particle_dynamics`          | `neighbor_search`                                                                                                                              | `all_pairs`, `cell_list`, or `auto`. The deterministic reference path is `O(N^2)` all-pairs. `cell_list` uses a sparse grid with a cell width covering the maximum pair interaction distance. `auto` retains all-pairs below 1024 active particles and selects the cell list at or above that measured crossover guard. |
| `particle_dynamics`          | `store_full_diagnostics`                                                                                                                       | Saves fluid velocity, angular velocity, rotation angle, collision force/count, gap ratio, wall weight, and two-wall warning arrays.                                                                                                                                                                                                 |
| `particle_dynamics`          | `two_wall_warning_gap_ratio`                                                                                                                   | Diagnostic threshold for flagging states where the opposite-wall gap is also small relative to the bubble radius. It does not change the single-wall model.                                                                                                                                                                         |
| `molecular_target_selection` | `default_mode`                                                                                                                                 | `manual` leaves the final choice to the Trame selector. `automatic` creates a deterministic synthetic target during candidate preparation; it does not infer a tumour.                                                                                                                                                          |
| `molecular_target_selection` | `influence_region_endothelial_wall_area_fraction`                                                                                              | Required in automatic mode and strictly in`(0,1)`. It sets the desired accessible wall area inside the compact physical influence region as a fraction of the complete theoretical network endothelial area.                                                                                                                      |
| `molecular_target_selection` | `target_positive_wall_fraction_within_influence`                                                                                               | Required in automatic mode and in`(0,1]`. It is the physical-area-weighted fraction of accessible wall inside the influence region that becomes target-positive.                                                                                                                                                                  |
| `molecular_target_selection` | `target_correlation_length_um`                                                                                                                 | Required positive physical correlation length of the patchy random field. It must span at least one CFD grid cell; values below three cells also trigger a resolution warning.                                                                                                                                                      |
| `molecular_target_selection` | `random_seed`, `random_field_modes`                                                                                                          | Reproducibility seed and Random Fourier Feature count. The seed is a non-negative integer and the mode count is at least 64; the shipped count is 512.                                                                                                                                                                              |
| `molecular_target`           | `enabled`, `region_mode`                                                                                                                     | Enables a target-positive subset of the true solid wall.`region_mode` is `none` or `mask_npz`; an enabled target cannot use `none`. The removed `physical_polygon` mode is rejected.                                                                                                                                      |
| `molecular_target`           | `mask_npz_path`, `mask_array_key`, `x_coordinates_key`, `z_coordinates_key`                                                              | Coordinate-aware Boolean ROI source for`mask_npz`. The named X/Z arrays are physical micrometer axes and must be one-dimensional, finite, and strictly increasing.                                                                                                                                                                |
| `molecular_target`           | `target_density_molecules_per_m2`                                                                                                              | Fixed target-site density on selected wall samples in`molecule/m^2`. A zero value is valid for a geometry-only no-bond contact pilot; enabled binding requires a positive value.                                                                                                                                                  |
| `molecular_binding`          | `enabled`, `model`                                                                                                                           | Enables the Revised-v7 deterministic mean-field bond extension. The only accepted model is`deterministic_mean_field`; it does not create discrete random bond identities.                                                                                                                                                         |
| `molecular_binding`          | `ligand_density_molecules_per_m2`                                                                                                              | Fixed ligand density on the bubble surface in`molecule/m^2`.                                                                                                                                                                                                                                                                      |
| `molecular_binding`          | `capture_distance_um`, `rest_length_um`                                                                                                      | Physical capture range and unstretched bond length in micrometers. These are separate geometric quantities and must be positive when binding is enabled.                                                                                                                                                                            |
| `molecular_binding`          | `association_rate_m2_per_molecule_s`                                                                                                           | Effective two-surface association coefficient in`m^2/(molecule s)`. A literature forward rate reported only in `s^-1` is dimensionally incompatible and must not be copied into this field.                                                                                                                                     |
| `molecular_binding`          | `zero_force_dissociation_rate_s`                                                                                                               | Zero-force first-order dissociation rate in`s^-1`; the Bell load factor modifies it at runtime.                                                                                                                                                                                                                                   |
| `molecular_binding`          | `bond_stiffness_pn_per_um`, `reactive_compliance_nm`, `temperature_k`                                                                      | Hookean bond stiffness, Bell reactive compliance, and absolute temperature. Their units are respectively`pN/um`, `nm`, and `K`.                                                                                                                                                                                               |
| `molecular_binding`          | `mean_field_warning_count`, `bell_exponent_limit`                                                                                            | Diagnostic threshold for low expected ligand/target counts and numerical cap for the Bell exponent. The count threshold is a warning, not a physical binding cutoff.                                                                                                                                                                |
| `binding_scenario_sweep`     | `enabled`                                                                                                                                      | Enables an empty-start no-bond exposure analysis and requires`molecular_binding.enabled=false`. The scenario table is built from pre-run inputs and is independent of whether this particular trajectory contains target contact.                                                                                                 |
| `binding_scenario_sweep`     | `da_on_reference_time_s`                                                                                                                       | Required finite positive reference time when the sweep is enabled. It defines the declared dimensionless scenario through`k_on = Da_on / (rho_T * t_ref)` and must be fixed before transport; it is never replaced by observed pilot contact time.                                                                                |
| `binding_scenario_sweep`     | `da_on_levels`, `capture_distance_to_rest_length_ratios`                                                                                     | Dimensionless on-rate sensitivity levels and reaction-geometry ratios used to create scenario reports; they are not fitted parameters.                                                                                                                                                                                              |
| `binding_scenario_sweep`     | `target_density_molecules_per_um2_levels`, `ligand_density_molecules_per_um2_levels`                                                         | Positive sensitivity axes expressed explicitly in`molecule/um^2`; the report converts them by `1 molecule/um^2 = 10^12 molecule/m^2`.                                                                                                                                                                                           |
| `cardiac_pulsatility`        | `enabled`                                                                                                                                      | Enables time-dependent carrier-flow and concentration-flux modulation. Missing sections default to disabled for backward-compatible programmatic tests; the production YAML enables it.                                                                                                                                             |
| `cardiac_pulsatility`        | `waveform`, `bpm`                                                                                                                            | Selects the positive`synthetic_ecg_envelope` surrogate and its cardiac period `T=60/bpm`. The waveform is not a measured pressure or velocity trace.                                                                                                                                                                            |
| `cardiac_pulsatility`        | `pulse_propagation_velocity_um_s`                                                                                                              | Converts root-to-grid vascular path distance to physical retarded time. The current production value is`25000 um/s`.                                                                                                                                                                                                              |
| `cardiac_pulsatility`        | `initial_phase_fraction`                                                                                                                       | Cardiac phase at formal empty-lumen time zero, expressed in`[0,1)` cycles.                                                                                                                                                                                                                                                        |
| `cardiac_pulsatility`        | `waveform_samples_per_cycle`                                                                                                                   | Number of samples in the compact single-cycle lookup table. Periodic linear interpolation is used at arbitrary internal-stage times.                                                                                                                                                                                                |
| `cardiac_pulsatility`        | `preserve_cycle_mean_flow`                                                                                                                     | Divides the positive surrogate by its cycle mean so the accepted steady CFD field and steady inlet rate remain cycle-mean references.                                                                                                                                                                                               |
| `cardiac_pulsatility`        | `modulation_strength`                                                                                                                          | Dimensionless pulse-amplitude control in`[0,1]`, applied as `a_eff = 1 + modulation_strength * (a - 1)` after optional mean normalization. `0` gives steady flow and inlet flux; `1` retains the full waveform. The shipped YAML uses `0.70`; `1.0` is only the loader fallback when this field is omitted.             |
| `red_blood_cell_transport`   | `root_discharge_hematocrit`                                                                                                                    | Root-vessel discharge haematocrit for the RBC-induced drift–diffusion reduced-order transport model. The formal default is `0.35`; `0.0` strictly disables drift, Fick correction, and stochastic dispersion. The RBC major diameter used only for the scale activation is fixed at `8 um` and is not configurable.                                                                        |
| `quick_test`                 | `enabled`                                                                                                                                      | Enables the same in-memory quick-test override path as the`--quick-test` command-line flag. Current default is `false`.                                                                                                                                                                                                         |
| `quick_test`                 | `grid_spacing_um`, `max_iterations`, `min_iterations`, `pressure_solver_max_iterations`, `n_steps`, `dt_s`, `integration_substeps` | Overrides applied when quick-test mode is enabled.                                                                                                                                                                                                                                                                                  |

The current configuration loader rejects legacy particle-contact controls with
a migration error. They are not silently mapped to the two current controls,
because a displacement-clipping parameter has different physical meaning from
predictive endpoint residual handling or true physical-time refinement.

### 4.1 Shipped YAML Effective Profile

The current checked-in YAML resolves to the following main-run profile before
any geometry-dependent inlet rate or population is known:

| Quantity                         | Current effective value                                                                                                           |
| -------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Grid spacing                     | `1.5 um`                                                                                                                        |
| Inlet concentration              | `2.0e7 MB/mL`                                                                                                                   |
| Formal output                    | `1000` intervals at `0.001 s`, hence `1001` frames over `1.0 s`                                                           |
| Particle integration             | `80` nominal subintervals per output interval, `1.25e-5 s` internal step, `80000` nominal subintervals total                |
| Revised-v15 contact controls     | `0.0015 um` endpoint residual-projection tolerance and at most `12` physical-time refinement levels                           |
| Diameter prior                   | Uniform over`1.5..2.5 um`; the scheduled marginal is subsequently weighted by finite-size inlet flux                            |
| Cardiac surrogate                | Enabled,`300 bpm`, `0.2 s` period, five formal cycles, mean preservation enabled, strength `0.70`                           |
| Current cardiac multiplier range | Approximately`0.71538..1.54581` for the configured 2048-sample cycle; periodic mean is one to floating-point precision          |
| Molecular extension              | Target, binding, and no-experiment sweep all disabled                                                                             |
| Candidate target selection       | Manual by default; automatic influence size, positive fraction, and correlation length are all unset rather than silently assumed |

The current quick profile has `41` stored frames over `0.04 s`, coarse
`8.0 um` spacing, two particle substeps per output interval, and
`internal_dt_s=5e-4 s`. It leaves concentration, diameter, cardiac, molecular,
and output sections unchanged.

Some important checked-in values differ intentionally from the fallback used
when a field is absent:

| Field                                       | Shipped YAML | Loader fallback if omitted |
| ------------------------------------------- | ------------ | -------------------------- |
| `domain.min_resolved_diameter_cells`      | `10.0`     | `8.0`                    |
| `field.flux_tolerance`                    | `1e-3`     | `1e-4`                   |
| `particle_dynamics.integration_substeps`  | `80`       | `1`                      |
| `cardiac_pulsatility.enabled`             | `true`     | `false`                  |
| `cardiac_pulsatility.modulation_strength` | `0.70`     | `1.0`                    |

### 4.2 Quick-Test, Path, and Timestamp Semantics

Quick mode is the logical OR of CLI `--quick-test` and YAML
`quick_test.enabled`. There is no CLI switch that disables a YAML-enabled quick
mode. The flag does not contain a second hard-coded config; it activates only
override keys actually present in the YAML `quick_test` mapping, and only the
seven fields listed in the configuration table above can be overridden.

Overrides modify temporary section mappings, not `cfg.raw`. Consequently a
CLI-activated quick run can save `run_config.yaml` with the original full-run
values and `quick_test.enabled: false`. Effective grid/time values are partly
available in `domain_metadata.yaml`, but the exact command line must also be
archived.

| Path-like input                    | Base used for a relative value                           |
| ---------------------------------- | -------------------------------------------------------- |
| CLI`--config`                    | Current working directory, followed by`Path.resolve()` |
| `input.swc_dir`                  | Repository root                                          |
| `output.results_dir`             | `ulm_microbubble_traj_gen` package directory           |
| `molecular_target.mask_npz_path` | Directory containing the selected YAML file              |

The timestamp is created during `load_config()`, not at first file write. The
default format has one-second resolution, and the runner creates the directory
with `exist_ok=True`. Two configs loaded in the same second can therefore point
to the same directory and overwrite like-named files or retain stale extras.
Serialize launches or use a more discriminating `timestamp_format`.

### 4.3 Loader Validation Boundaries

- Use native YAML booleans. The loader calls Python `bool(...)`; a quoted value
  such as `"false"` is a non-empty string and is therefore truthy.
- Particle dynamics always validate integrator, substep count, mobility blend,
  collision, and search settings. Cardiac waveform name, BPM, propagation
  speed, phase, sample count, and modulation strength are validated even when
  cardiac pulsatility is disabled.
- Molecular target mode, density, and mask key names are checked at config
  load. A mask path is required for enabled `mask_npz` mode, while filesystem
  existence, arrays, axes, and shape are validated during target construction.
  An ROI selecting zero eligible wall
  sites is allowed, so inspect the saved target count before interpreting a run.
- Molecular target selection mode is `manual` or `automatic`. Automatic mode
  requires an influence-region wall-area fraction in `(0,1)`, a positive-wall
  fraction in `(0,1]`, a positive physical correlation length, a non-negative
  integer seed, and at least 64 random-field modes. Revised-v9 keys
  `target_endothelial_wall_area_fraction` and `maximum_disconnected_regions`
  are rejected with an explicit migration error instead of being reinterpreted.
- Enabled binding requires an enabled target with positive target density and
  positive ligand density, capture/rest lengths, effective 2D association
  coefficient, stiffness, and reactive compliance. The zero-force off-rate may
  be zero. Ambiguous paper-style one-dimensional rate/encounter-radius aliases
  are rejected even when binding is disabled.
- An enabled sensitivity sweep requires target enabled, binding disabled,
  positive rest length, full diagnostics, contact pilot and observed-contact
  guard enabled, and four non-empty positive level groups.
- Config loading directly validates only positive finite particle
  concentration. Diameter bounds, geometry/path controls, capacity, inlet
  feasibility, and capture distance versus each sampled bubble are checked in
  downstream builders. Input file existence is also deferred to the runner.

### 4.4 Grid Resolution Note

`grid_spacing_um` controls mask quality, flow cost, velocity-gradient accuracy,
continuous wall-distance accuracy, and therefore near-wall mobility. The current
full-run default is `1.5 um`; the current quick-test default is `8.0 um`.

The formal narrow-cell criterion is:

```text
D_eff_px < 8
```

The generation-stage warning threshold is stricter:

```text
2 * min(vessel_radius_um) / grid_spacing_um < 10
```

This warning does not stop the run by itself. It indicates that the smallest
physical vessels lack the recommended rasterization margin, which is especially
important for wall-shear- and mobility-sensitive runs. In the current default
model, the smallest exported radius is about `7.642 um`; its diameter spans
about `10.19` cells at `1.5 um`, just above the configured `10`-cell warning
threshold. The `8.0 um` quick test resolves the same diameter with only about
`1.91` cells and is therefore a plumbing test, not a wall-sensitive production
setting.

## 5. Lumen Mask Generation

The active rasterizer follows the v4 geometry-first policy: preserve DCCO
topology and physical radii, repair geometry before rasterization, and avoid
global dilation/closing as a formal fix.

### 5.1 Continuous Geometry

Each vessel segment is projected into the X-Z plane:

```text
p0 = (x_p, z_p)
p1 = (x_d, z_d)
segment = p1 - p0
length = |segment|
```

Each segment becomes a corridor polygon generated from its centerline and a
geometric radius:

```text
mask_radius = max(physical_radius, min_lumen_radius_cells * grid_spacing_um)
```

For the default model and current default spacing, this fallback does not
inflate the smallest vessels because the physical radii are larger than the
minimum mask radius. In general, the fallback should be treated as a raster
visibility guard; physical diagnostics still use the exported DCCO radius.

Open endpoints are treated differently from internal junctions:

```text
root inlet and terminal outlet endpoints: flat cap
internal junction endpoints: explicit junction geometry
```

### 5.2 Junction Geometry

Internal junctions are reinforced in continuous geometry before rasterization.
For each non-open endpoint shared by multiple incident vessels, the code builds:

```text
junction disk
+ short overlap corridors along incident vessel directions
+ a smooth convex-hull junction polygon
```

The junction core radius is:

```text
R_j = R_max + 2 * grid_spacing_um
```

The incident overlap length is:

```text
L_j = max(2 * R_j, 6 * grid_spacing_um)
```

This reduces single-pixel holes and narrow necks at bifurcations without
changing the underlying vessel graph or inventing new vessel segments.

### 5.3 Polygon Union and Hole Removal

All corridor and junction polygons are merged with a geometric union. The
merged polygon is repaired if it is invalid, and interior polygon rings are
removed before rasterization. This makes `holes=0` primarily a continuous
geometry property instead of a late image-processing artifact.

### 5.4 Exact Shapely Cell Coverage

For every Cartesian cell \(C_{ij}\), the rasterizer evaluates the exact
intersection area with the continuous lumen \(\Omega\):

\[
\alpha_{ij}
=\frac{\operatorname{Area}(C_{ij}\cap\Omega)}
       {\operatorname{Area}(C_{ij})}.
\]

The saved `lumen_fraction` stores \(\alpha_{ij}\). A cell enters the Boolean
sampling mask when \(\alpha_{ij}\ge0.5\) and its centre lies inside the
continuous lumen. The centre check guarantees that every accepted DOLFINx
sampling point belongs to the boundary-fitted finite-element mesh.

After rasterization, only conservative cleanup is allowed:

```text
fill internal binary holes
remove isolated fragments of at most 4 cells
```

No global dilation, global closing, or erosion is used as a formal repair step.

### 5.5 Segment Attribute Assignment

After the lumen mask is fixed, per-cell vessel attributes are assigned by the
nearest centerline among physical vessel segments. This stage does not expand
the mask. It fills:

```text
vessel_id
radius_um
flow_rate_um3_s
q2d_flow_um2_s
viscosity_mpas
direction_xz
distance_to_centerline_um
```

Junction-added cells that are not directly assigned by a centerline inherit
attributes from the nearest assigned lumen cell. The separate
`junction_core_mask` records the junction geometry so narrow diagnostics can
avoid treating a junction core as if it belonged exclusively to one narrow
child branch.

### 5.6 Wall Fields

Grid wall fields are sampled from the same continuous solid wall used by
DOLFINx and particle contact:

```text
d_ij, n_ij = continuous_geometry.exact_solid_wall_state(x_ij)
wall_mask = lumen_mask AND d_ij <= grid_spacing_um / sqrt(2)
```

The saved `distance_to_wall_um` is \(d_{ij}\), and
`wall_normal_xz` is the continuous inward normal \(\mathbf n_{ij}\). No binary
distance transform or pixel-gradient normal reconstruction is used.

### 5.7 Hybrid Particle Hydrodynamics and Continuous Wall Geometry

Particle transport uses one hybrid representation. Velocity, velocity
gradient, viscosity, local radius, and wall-shear proxy are interpolated from
the regular CFD grid. Wall distance, inward normal, contact, swept-path
acceptance, inlet/outlet crossing, and finite-radius accessibility are queried
only from the required pre-raster `ContinuousVesselGeometry`.

Let \(\Omega\) be the continuous lumen, let \(\Gamma_o\) be its anatomical
inlet and outlet sections, and let the solid wall be

\[
\Gamma_w = \partial\Omega\setminus\Gamma_o,
\qquad
\Gamma_w\cap\Gamma_o=\varnothing.
\]

The continuous object supplies lumen membership, wall distance and inward
normal, directed open sections, target-wall support, and finite-radius inlet
connectivity. Particle transport does not reconstruct these decisions from
`wall_mask`, raster labels, a grid distance transform, or a grid normal field.
The sampled Boolean wall/open overlays are retained only for indexing,
selection, and visualization; they do not determine particle mechanics.

For radius \(R\), the true finite-size gap is

\[
\boxed{g_R(\mathbf x)=d_w(\mathbf x)-R}.
\]

This continuous definition is shared by inlet admission, swept-path
acceptance, predictive contact, saved `record_wall_gap_um`, and molecular
capture. Every accepted continuing centre must satisfy \(g_R\ge0\). The
configured contact tolerance limits only a small endpoint residual projection;
it does not permit overlap and is not a physical stand-off.

The hydrodynamic mobility uses a separately regularized gap only when its
near-wall coefficients are evaluated:

\[
h_{\mathrm{hyd}}=\max(g_R,h_{\mathrm{reg}}).
\]

In the implementation, `xi_min` provides this algebraic lower scale through
the dimensionless gap used by the mobility formula. It never changes
\(g_R\), \(\Omega\), target reaction geometry, or lifecycle acceptance. This
separation prevents the mobility regularizer from becoming an artificial
positive-clearance wall.

For each radius, the continuous geometry evaluates the inlet-connected
accessible component

\[
\Omega_R^{\mathrm{in}}
=\operatorname{Conn}_{\mathrm{in}}
\{\mathbf x\in\Omega:d_w(\mathbf x)\ge R\}.
\]

This prevents a finite-size bubble from entering an inlet-disconnected or
too-narrow branch without introducing a second raster accessibility model.
For a complete straight trial chord, the swept disc is checked against the
continuous boundary representation. A forbidden pocket between two feasible
endpoints therefore requests
real-time refinement instead of being hidden by endpoint-only testing.

## 6. Mask Quality and Narrow-Cell Diagnostics

The run stops before flow solving unless the rasterized lumen is one connected
fluid domain and has no internal holes:

```text
components = 1
holes = 0
```

Narrow cells are not near-wall cells. They are lumen cells whose effective
local diameter is under-resolved:

```text
D_graph_px = 2 * radius_um / grid_spacing_um
D_mask_px  = 2 * medial_axis_distance_px propagated to lumen cells
D_eff_px   = min(D_graph_px, D_mask_px)
narrow     = D_eff_px < 8
```

Special diagnostic rules:

```text
junction_core cells: use D_mask_px for narrow classification
open-boundary diagnostic region: use D_graph_px to avoid treating an open cap as a closed wall
```

The diagnostics report:

```text
narrow_cells
narrow_fraction
narrow_in_inlet_or_outlet
narrow_in_junction_core
narrow_in_main_junction_core
junction_core_cells
min_resolved_diameter_px
p1_resolved_diameter_px
p5_resolved_diameter_px
median_resolved_diameter_px
```

`narrow_segments.csv` additionally records segment-level radius, graph
diameter, mask diameter, effective diameter, length, flow, terminal status, and
whether the segment is an inlet/outlet segment.

## 7. Open-Edge Boundary Conditions

The solver uses mask-boundary faces as open boundaries. It does not use thick
interior velocity bands for inlet or outlet conditions.

For each root vessel, inlet faces are selected near the proximal endpoint. For
each terminal vessel, outlet faces are selected near the distal endpoint. Faces
are selected from lumen cells adjacent to non-lumen space and must be aligned
with the intended outward normal. Each selected face has:

```text
cell index
face index
axis
normal
face length
label
kind: inlet or outlet
signed flux
```

The 3D exported flow rate is converted to a 2D flux target by:

```text
q2d_um2_s = flow_rate_um3_s / effective_thickness_um
```

Within one inlet or outlet, flux is distributed over selected faces with a
parabolic lateral weight multiplied by face-normal alignment. Final actual
fluxes are always integrated from the final face flux arrays, not from a
separate diagnostic approximation.

The accepted `FlowField` carries these same selected faces as authoritative
`open_face_*` arrays and groups them into oriented `open_section_*` records.
Particle geometry copies those arrays, removes exactly those faces from
\(\Gamma_w\), and uses the outlet section point, outward normal, tangent, and
half-width to evaluate \(\psi_k\). This is the required identity between CFD
open boundaries and particle lifecycle geometry; cell labels are not rebuilt
into a second set of openings on the production path.

## 8. Velocity Solver

The active velocity solver is a hybrid method:

```text
PhiFlow explicit diffusion for viscous relaxation
SciPy sparse CG for finite-volume pressure projection
```

It is intentionally not labeled as a pure PhiFlow incompressible solver. The
canonical solver identity is:

```text
phiflow_viscous_fv_projection_2d
```

The solver rejects the misleading name:

```text
phiflow_incompressible_2d
```

### 8.1 Initial Velocity

The initial velocity field uses per-cell 2D flux, radius, centerline distance,
and a smoothed local vessel direction. The profile is Poiseuille-like in the
2D lumen mask:

```text
relative = distance_to_centerline_um / radius_um
profile  = max(0, 1 - relative^2)
speed    = 3 * q2d / (4 * radius_um) * profile
velocity = smoothed_direction_xz * speed
```

### 8.2 Iteration Loop

Each solver iteration performs:

```text
1. convert current centered velocity to a PhiFlow staggered grid
2. apply explicit viscous diffusion when viscosity and dt are positive
3. convert the diffused centered velocity to face fluxes
4. enforce fixed open-edge inlet/outlet face fluxes
5. project internal fluid-fluid face fluxes with a pressure solve
6. reconstruct the cell-centered velocity from projected face fluxes
7. compute physical convergence metrics from the same final face fluxes
```

A tqdm progress bar prints during this loop:

```text
Solving 2D lumen flow
```

### 8.3 Finite-Volume Projection

Each lumen cell is a fluid control volume. Internal fluid-fluid faces are used
to build the sparse pressure matrix. Solid-wall faces have fixed zero normal
flux. Open-boundary faces have fixed prescribed inlet or outlet flux.

The per-cell net outflow is:

```text
net = flux_x[i+1,j] - flux_x[i,j] + flux_z[i,j+1] - flux_z[i,j]
```

The finite-volume divergence is:

```text
divergence_s_inv = net / grid_spacing_um^2
```

The pressure correction is solved by SciPy CG, with an optional PyAMG smoothed
aggregation preconditioner when available. The exported cell-centered velocity
is reconstructed from the final projected face fluxes:

```text
u_x = 0.5 * (flux_x_left + flux_x_right) / grid_spacing_um
u_z = 0.5 * (flux_z_down + flux_z_up) / grid_spacing_um
```

The saved velocity stage is:

```text
face_flux_after_projection
```

No wall cleanup or inlet/outlet overwrite is allowed after the final projection.

### 8.4 Physical Acceptance Criteria

The steady momentum-balance check uses the actual velocity increments produced
by the implemented split scheme. Let `u^k` be the previous projected velocity,
`u^nu` the canonical face-flux-reconstructed velocity after viscous relaxation,
`u^b` the velocity after applying the fixed open-boundary fluxes, and
`u^(k+1)` the final projected velocity. On the interior lumen mask, which
excludes the open-boundary band, define:

```text
a_viscous  = (u^nu - u^k) / viscous_relaxation_dt_s
a_pressure = (u^(k+1) - u^b) / viscous_relaxation_dt_s
r_momentum = a_viscous + a_pressure

normalized_momentum_residual =
    norm(r_momentum) / (norm(a_viscous) + norm(a_pressure) + epsilon)
```

This is the residual of the implemented discrete steady Stokes-type balance.
It does not reinterpret the saved projection potential as physical pressure,
and it does not add an inertial or convective Navier-Stokes term. The complete
lumen update residual, including the prescribed open-boundary increment, is
saved separately for diagnostics and is not the interior acceptance quantity.

The solver accepts a velocity field only if all conditions hold after the
minimum iteration count:

```text
pressure_cg_info == 0
final_normalized_divergence_error <= divergence_tolerance
interior_normalized_divergence_error <= divergence_tolerance
max(final total/target/per-label flux errors) <= flux_tolerance
wall_penetration_max_um_s <= max_wall_penetration_um_s
normalized_momentum_residual is finite
normalized_momentum_residual <= momentum_residual_tolerance
```

If these criteria fail, `FlowConvergenceError` is raised. The failed field and
diagnostics may be written, but particle advection is stopped.

Accepted reusable fields must also carry the
`discrete_stokes_momentum_balance_v1` acceptance schema, a true
`momentum_residual_converged` flag, and a finite residual satisfying the current
configuration. Historical fields accepted before this schema cannot bypass the
new momentum gate.

The checked-in production configurations allow up to `200` outer iterations,
and quick-test mode allows up to `120`. These are safety ceilings rather than
fixed work counts; the solve still stops immediately after every acceptance
condition passes. The larger ceilings prevent the momentum gate from being
defeated by the former mass-conservation-only early exit.

## 9. Wall-Shear Proxy

The current WSS field is an in-plane strain-rate proxy, not a full 3D CFD wall
shear reconstruction. It uses the final accepted velocity, per-cell viscosity,
and wall-normal/tangent vectors derived from the mask distance field.

Velocity gradients are computed on the X-Z grid. The code forms the symmetric
strain contribution along the local wall normal and projects it onto the local
tangent. The magnitude is multiplied by apparent viscosity:

```text
shear_rate = tangent . (strain_tensor * normal)
wall_shear_stress_pa = abs(viscosity_pa_s * shear_rate)
```

Outside-lumen WSS is set to zero.

The mobility model does **not** use this nonnegative WSS proxy as its driving
shear. It bilinearly samples the complete precomputed velocity-gradient tensor
and retains the signed local shear rate

```text
S = tangent.T @ grad(velocity) @ inward_normal
```

so tangential translation and rotation retain their physical direction. The
existing WSS proxy is still sampled and saved for trajectory diagnostics and
visualization.

## 10. Particle Initialization and Advection

### 10.1 Model Dispatch and Convergence Guard

`runner.py` accepts the flow result, precomputes mobility and cardiac fields,
and directly invokes the only supported continuous-perfusion particle model.
The solver must report:

```text
flow_field.solver_metadata["physical_converged"] == true
```

to be true before the runner reaches particle transport. Particle transport
never runs on an unaccepted field. The former particle-mode selector and
`passive_tracer` path have been removed.

```text
accepted flow
    -> precompute hydrodynamic fields
    -> construct continuous-perfusion inlet flux
    -> save the empty lumen as formal frame 0
    -> synchronous finite-size mobility transport
```

The fixed-population, passive-tracer, and Taichi particle implementations have
been removed. Positive inlet concentration is validated when configuration is
loaded, and the generator has a single continuous-perfusion call chain.

### 10.2 Deterministic Number Flux, Capacity, and Permanent IDs

For inlet coordinate `s`, radius `R`, inward section normal `n_in`, effective
thickness `b`, and physical concentration `C0`, the implemented joint flux is:

```text
u_n_plus(s) = max(velocity(s) dot n_in, 0)
j(s,R) = C0 * b * u_n_plus(s) * p_R(R) * I[s in Omega_R^in]
lambda = integral_R integral_s j(s,R) ds dR
```

The configured concentration is in `MB/mL`; multiplying by `1e6` converts it to
`bubble/m^3`. Length, velocity, section width, and effective thickness are in
micrometre units, so the section volume flux is multiplied by `1e-18` before it
is combined with `C0`. Between configured radius bounds, `p_R` is a uniform
prior. The inverse marginal CDF is nevertheless weighted by finite-size
accessible flux, so the scheduled radius population need not be uniformly
distributed when larger bubbles lose usable inlet width or encounter a
different velocity profile.

```text
t_k = (k + 1/2) / lambda
eta_k = Halton_base2(k + 1)
xi_k  = Halton_base3(k + 1)
R_k   = inverse marginal flux CDF(eta_k)
s_k   = inverse conditional position flux CDF(xi_k | R_k)
```

This sequence has no random seed. The midpoint event count obeys
`abs(N_in(T) - lambda*T) <= 1/2`, and changing output `dt_s` does not change the
physical event table. Diameter remains immutable for permanent ID `k`.

The schedule horizon is exactly the formal recording duration.
Its exact event count defines automatic state capacity when
`max_unique_bubbles=0`; an explicit smaller cap fails before transport. Formal
records are ragged, and the runtime checks their cumulative count against
`max_particle_frame_records`. At the end the log reports how many formal
frames that limit supports at the observed peak active population.

During formal trajectory generation, a live `tqdm` bar named
`Particle transport` advances once for every completed
**nominal** internal subinterval. Its total is exactly
`n_steps * integration_substeps`, so the percentage, processing rate, and ETA
remain consistent with the selected base temporal resolution. An injection
event can split one nominal subinterval into multiple accepted segments, but
those event-split segments do not increase the progress-bar total. At each
saved-frame boundary, the postfix reports the completed output frame, active
bubbles, and waiting events. There is no unrecorded equilibration bar. If
`tqdm` is unavailable, progress display silently degrades without changing
numerics.

Every physical microbubble receives a monotonically increasing permanent ID.
Its state array index equals that ID and is never reused. Formal frames contain
only current active IDs; no fixed-width inactive lanes are written.

### 10.3 Fixed Injection Section, Waiting, and Empty Formal Start

For every root vessel, one physical line section is placed downstream of the
open-boundary correction depth but within the straight root segment, before its
distal junction. The accepted velocity and solid-wall distance are sampled on
that fixed section. The code logs the integrated section flow and its difference
from the accepted inlet/reference flow.

Every scheduled centre sample is required to lie inside the canonical lumen,
have \(g_R\ge0\), and belong to \(\Omega_R^{\mathrm{in}}\) for its assigned
radius. Radius-dependent accessible section flux therefore controls both the
total number flux and the sampled diameter/position distribution. A candidate
that violates any canonical geometry condition is rejected rather than moved
to a nearby raster cell.

At `t_k`, the scheduled bubble keeps its unique `(R_k, s_k)`. It enters only if:

```text
center_distance >= R_k + R_j  for every active j
```

Otherwise it remains in a deterministic upstream queue. It is not moved,
resized, resampled, deleted, or force-inserted. Waiting entries are retried after
accepted integration segments; their actual admission and waiting times are
stored separately from planned time. Retry scans preserve deterministic ID
order, but admission is **not strict FIFO**: a later ID may enter at its own
non-overlapping inlet location while an earlier ID remains blocked. Remaining
blocked entries keep their relative order for the next scan.

Formal frame 0 is the empty lumen at physical `t=0`; there is no unrecorded
warm-up and no inlet/outlet or population-equilibrium gate. The first possible
injection remains the midpoint event `t_0=1/(2*lambda)`. Consequently, the
saved trajectory intentionally includes the physical filling transient and its
age-coherent transport front. If the formal horizon is shorter than `t_0`, the
writer produces a valid sequence of empty frames rather than failing.

### 10.4 Hydrodynamic Field Precomputation

The mobility path builds the following grid fields once per accepted flow:

```text
velocity_gradient_s_inv[..., velocity_component, coordinate]
    = [[du_x/dx, du_x/dz],
       [du_z/dx, du_z/dz]]

vorticity_y_s_inv = du_x/dz - du_z/dx
dynamic_viscosity_pa_s
solid_wall_distance_um
inward_wall_normal_xz
solid_site_mask
open_boundary_mask
canonical boundary geometry: Omega, Gamma_w, open faces, directed sections psi_k
inlet_connected_clearance_um for all radius queries
```

Because velocity is `um/s` and coordinates are `um`, every gradient component
has units `1/s`. Raster viscosity is converted from `mPa s` to `Pa s`. Invalid
or zero-valued solid viscosity cells are filled from the nearest positive lumen
cell before interpolation, preventing a near-wall 2x2 stencil from being
diluted by the solid storage value.

All particle positions, matrix coefficients, pair forces, and integration
states use `float64` internally. Grid storage and the pre-existing trajectory
fields may remain `float32`; the molecular state, rate, force, torque, and
reaction-area arrays are deliberately written as `float64` so a finite
double-precision Bell/formation diagnostic cannot overflow merely because it
was down-cast during output.

### 10.5 Local Coordinates and Unit Convention

The mobility implementation retains planar translation and out-of-plane
rotation:

```text
U = (V_x, V_z, Omega_y)
G = (F_x, F_z, T_y)
```

For inward unit normal `n = (n_x, n_z)`, the tangent is fixed as:

```text
t = (-n_z, n_x)
t cross n = +Y
```

The scaled local variables are:

```text
U_tilde = (V_t, V_n, R * Omega_y)
G_tilde = (F_t, F_n, T_y / R)
```

The numerical unit system is:

```text
length       um
time         s
velocity     um/s
viscosity    Pa s
force        pN
torque       pN um
angular rate rad/s
```

With these units, `1 / (6*pi*mu*R)` numerically maps `pN` to `um/s` without
repeated powers-of-ten conversions inside the high-frequency kernels.

### 10.6 Free-Space and Near-Wall Mobility

The free-space scaled mobility is:

```text
M_hat_infinity = diag(1, 1, 3/4)
```

The true geometric gap and hydrodynamic evaluation gap are deliberately
distinct:

```text
g_R     = solid_wall_distance - R
h_hyd   = max(g_R, 0)
xi      = h_hyd / R
xi_eff  = max(xi, xi_min)
```

`g_R` is the canonical wall gap used by geometry and is saved directly. A
negative accepted value is a Revised-v15 invariant failure. `h_hyd` and
`xi_eff` are algebraic guards used only by the mobility coefficients; neither
is a permitted centre clearance or a molecular capture gap. The
supplementary-material resistance coefficients are implemented as:

```text
C_Ft = (8/15) * ln(xi_eff) - 0.9588
C_Tt = (1/10) * ln(xi_eff) - 0.1895
C_Fr = (2/15) * ln(xi_eff) - 0.2526
C_Tr = (2/5)  * ln(xi_eff) - 0.3187

Delta  = C_Tt*C_Fr - C_Ft*C_Tr
Lambda = 1 + 1/xi_eff
```

The scaled single-planar-wall mobility is:

```text
M_hat_wall =
    [[C_Tr/Delta, 0, 3*C_Fr/(4*Delta)],
     [0, 1/Lambda, 0],
     [C_Tt/Delta, 0, 3*C_Ft/(4*Delta)]]
```

The two published translation-rotation cross terms are retained separately;
the implementation does not silently symmetrize them. Their maximum relative
reciprocity mismatch is saved as a diagnostic.

The near-wall asymptote is not extrapolated through the vessel center. Its
evaluation gap is frozen at `xi_near`, then blended into free space:

```text
xi_star = min(max(xi, xi_min), xi_near)
s       = clip((xi - xi_near) / (xi_far - xi_near), 0, 1)
w       = 1 - 3*s^2 + 2*s^3

M_hat_eff = (1 - w)*M_hat_infinity + w*M_hat_wall(xi_star)
```

The default `xi` values are engineering transition parameters and must not be
reported as universal physical constants.

### 10.7 Background Hydrodynamic Motion

The bulk background motion is:

```text
V_bulk     = bilinear(velocity, position)
Omega_bulk = 0.5 * (du_x/dz - du_z/dx)
```

Near a wall, the signed shear rate is:

```text
S = t.T @ grad(velocity) @ n
```

The supplementary-material equivalent shear load uses two distinct heights:

```text
h_gap    = h_hyd
H_center = R + h_gap

F_t_shear = 1.7005 * 6*pi*mu*R*H_center*S
T_y_shear = 0.9440 * 4*pi*mu*R^3*S
```

The force coefficient uses the Goldman, Cox, and Brenner (1967) Part II
contact-limit value `F_s = 1.7005`. The `1.0075` value printed in the downstream
Borden supplementary material is treated as a transcription error. The torque
coefficient remains the Goldman contact-limit value `T_s = 0.9440`.

Separating `h_gap` from `H_center` is essential: the shear force does not
incorrectly vanish when the bubble surface reaches the wall. The local
near-wall shear motion is obtained from `M_hat_wall`, and the complete
background hydrodynamic state is blended separately:

```text
U_H_tilde = (1 - w)*U_bulk_tilde + w*U_wall_shear_tilde
```

This avoids applying the transition weight twice.

#### 10.7A Cardiac carrier-flow modulation

The ECG-shaped helper is retained as a positive engineering surrogate rather
than being interpreted as an electrical-to-flow constitutive law. Only one
cycle is stored. When `preserve_cycle_mean_flow=true` (the shipped setting),
the sampled cycle is first normalized to mean one. Pulse strength then scales
only its deviation from the steady multiplier:

```text
T                  = 60 / bpm
a_base             = a_raw / mean_cycle(a_raw) if preserve_cycle_mean_flow=true
                     a_raw                       otherwise
a_eff              = 1 + modulation_strength*(a_base - 1)
tau(x)             = root_path_distance(x) / pulse_propagation_velocity
phase_time(x,t)    = t + initial_phase - tau(x)
a(x,t)             = periodic_linear_interpolation(a_eff, phase_time)
u(x,t)             = a(x,t)*u_steady(x)
```

If mean preservation is disabled, normalization is skipped and the resulting
cycle is not guaranteed to have unit mean. The shipped strength is `0.70`;
zero makes both carrier flow and inlet event rate steady, while one retains the
full surrogate deviation.

Root-path distance is obtained from vessel parent topology plus the projected
position within the attributed segment. It is a spatial fluid property, not a
particle's accumulated travel distance. Delay remains in physical seconds; no
conversion through samples-per-heartbeat is used. The delay and waveform are
sampled at every Euler stage and at both the current and predicted Heun stage.

For consistency with the retarded spatial multiplier, the velocity gradient
used by vorticity and near-wall mobility is reconstructed as:

```text
grad(u) = a*grad(u_steady) + u_steady tensor grad(a)
grad(a) = -a_time_derivative * grad(tau)
```

The saved nonnegative WSS proxy sampled along a trajectory is multiplied by the
same local positive factor. The multiplier therefore affects carrier velocity,
its reconstructed gradient/vorticity, near-wall shear motion, saved WSS proxy,
and the inlet number-flux clock. Collision force and molecular bond force are
not directly multiplied by `a`; they can still change indirectly because the
modulated flow changes relative particle geometry and bond slip. The current
implementation has no separate pulse-strength control for injection events.

The concentration-calibrated inlet event rate is time-dependent:

```text
lambda(t) = lambda_mean * a_inlet(t)
integral(0, t_k, lambda(t) dt) = k + 0.5
```

Event times are obtained by deterministic inversion of this cumulative number
flux. The radius/location Halton sequence, empty formal frame 0, deterministic
ordered-retry waiting queue, permanent IDs, and never-reused state slots are
unchanged. Over an integer number of mean-normalized cardiac cycles, the
integrated event count recovers the steady mean-rate count apart from the
existing half-event rule.

### 10.8 Pair Collision Model

For particles `i` and `j`:

```text
compression_ij = max(0, R_i + R_j + collision_layer_um - distance_ij)
n_ij           = (x_i - x_j) / distance_ij

mu_rel = n_ij.T @ (M_VF_i_global + M_VF_j_global) @ n_ij
k_ij   = 1 / (mu_rel * collision_relaxation_time_s)
F_ij   = k_ij * compression_ij * n_ij
F_ji   = -F_ij
```

Coincident centers use a deterministic direction derived from the two
permanent IDs. Relative mobility must be finite and strictly positive or the
run stops. Central collision force has no explicit torque, but a tangential
force can induce rotation through near-wall translation-rotation coupling.

All pair forces are evaluated from one immutable particle snapshot before any
position is committed. Therefore traversal order cannot make a later particle
react to an already-updated earlier particle. `all_pairs` is the deterministic
Numba-capable `O(N^2)` reference implementation. The sparse `cell_list` visits
only the current cell and its eight neighbours, orders equal-cell candidates by
permanent ID, and accumulates each particle's force without shared writes.
`auto` keeps all-pairs for fewer than 1024 active particles and uses the cell
list at larger populations. This threshold is a performance choice, not a
physical parameter. The public path has no fixed default population: cost must
be assessed against the observed peak active count produced by concentration,
inlet flux, residence time, exits, and waiting.

### 10.9 Unified Velocity, Contact Constraint, and Time Integration

Before the rigid-wall reaction is added, the generalized free velocity is:

```text
U0_tilde
    = U_H_tilde
    + (1 / (6*pi*mu*R)) * M_hat_eff
        * (G_collision_tilde + G_bond_tilde)
```

Collision torque is zero in the implemented central pair model; the optional
bond generalized load contains tangential/normal force and usually a nonzero
torque. The summed load passes through the local mobility operator once. Here
free means only that the rigid-wall reaction has not yet been added: carrier
flow, near-wall shear, pair collisions, and molecular loads are already in
`U0`.

Let \(\mathcal M\) be the complete physical 3-by-3 mobility mapping
\((F_x,F_z,T_y)\) to \((V_x,V_z,\Omega_y)\), let the local inward wall normal
be \(\mathbf n=(n_x,n_z)\), and define

\[
\mathbf a=(n_x,n_z,0)^T,
\qquad
\mathbf B=
\begin{pmatrix}1&0&0\\0&1&0\end{pmatrix}.
\]

For a complete positive physical interval \(\tau\), Revised v15 first predicts
the normal end-gap with the exact accepted start gap:

\[
\boxed{
g_{\mathrm{pred}}=g_R(\mathbf x_n)+
\tau\mathbf a^T\mathcal U^0
}
\]

and obtains the unique single-wall unilateral reaction in closed form:

\[
\boxed{
\lambda=\max\left(0,
-\frac{g_R(\mathbf x_n)/\tau+\mathbf a^T\mathcal U^0}
{\mathbf a^T\mathcal M\mathbf a}\right),
\qquad \mathbf a^T\mathcal M\mathbf a>0
}
\]

The constrained velocity and trial endpoint are then

\[
\mathcal U^c=\mathcal U^0+\mathcal M\mathbf a\lambda,
\qquad
\mathbf x^*=\mathbf x_n+\tau\mathbf B\mathcal U^c.
\]

The centered reaction has no direct torque component, but the full mobility
matrix may couple it into rotation. It adds no tangential friction and is not a
molecular bond force. There is no nonlinear root search, reaction-ray bracket,
PGS iteration, impulse, friction, or Baumgarte correction in the v15 path.

The authoritative distance \(d_w\), true gap, first-contact location, and
normal are analytic queries against the finite solid face segments in
\(\Gamma_w\); the raster node-distance field is retained for topology and
visualization only. The complete centre chord is checked as a swept disc. If a
trial endpoint has only a small residual \(-\epsilon\le g_R<0\), it may be
moved outward once along the exact distance gradient. A larger endpoint
penetration or any negative interior swept gap rejects the complete trial and
requests physical-time bisection. A simultaneous contact with distinct solid
faces is a hard model-limit error, not an invitation to choose or average one
normal.

`record_velocities_um_s` stores the last accepted internal v15 constrained
generalized translation for that permanent ID when one exists. A newly
admitted bubble without a prior accepted internal step uses the current
saved-state RHS as its fallback. It is neither plain fluid velocity nor the
saved-frame finite-difference velocity.

`record_realized_velocities_um_s` stores that same-permanent-ID accepted
position difference directly. The first observation uses the following
same-ID interval when available, interior observations use the following
interval, and the final observation reuses the preceding interval. Different
permanent IDs are never differenced even when a visualization lane is reused.

With full diagnostics, `record_fluid_velocities_um_s` stores the local carrier
velocity after any cardiac modulation but before particle wall, collision, or
bond mobility corrections. Without cardiac modulation it reduces to the direct
bilinear sample of the accepted steady CFD field.

`particles.dt_s` is the interval between saved observations. The mobility path
divides it into an integer number of exactly aligned internal steps:

```text
output_dt_s = particles.dt_s
internal_dt_s = output_dt_s / particle_dynamics.integration_substeps
```

Mobility, collisions, rotation, the predictive reaction, and optional molecular
loads are recomputed for every accepted internal interval. Only output-boundary
states are written to the trajectory file, so visualization cadence and file
size do not scale directly with the accepted substep count.

When a planned `t_k` lies inside an internal step, that step is split at the
event time. Neither resulting segment is longer than `internal_dt_s`; therefore
the planned event is not rounded to the output or internal **time** grid once
the inlet flux and waveform have been constructed. It is not independent of
the spatial CFD grid, because that grid affects the accepted inlet velocity and
finite-size accessibility. The selected Euler/Heun method still advances every
equal-time active-particle snapshot synchronously. A blocked scheduled position
is retried at accepted segment boundaries, so its admission delay remains a
substep-convergence quantity.

Euler evaluates one free RHS and passes it, the complete mobility, radius, and
complete interval to the v15 predictor. The predictor obtains its wall normal
from the exact first-contact point. Strict Heun first builds a feasible v15
predictor, recomputes carrier flow, mobility, collisions, and bonds there,
averages the two generalized velocities and mobilities, and solves the corrected
endpoint with the same closed-form v15 constraint. A contact normal is always a
geometry value and is never averaged as an ODE state. Speculative predictor bond state is discarded; molecular state is
committed only for the accepted corrected interval.

After either stage, the complete chord from \(\mathbf x_n\) to the proposed
endpoint must remain feasible. If an interior negative-gap pocket, invalid
normal mobility, or inconsistent stage makes it unacceptable,
the **whole physical interval** is replaced by

\[
[t_0,t_1]
\longrightarrow
[t_0,(t_0+t_1)/2]+[(t_0+t_1)/2,t_1].
\]

Both half intervals recompute the RHS, mobility, normal, reaction, outlet
events, rotation, and molecular state in chronological order. Their durations
sum exactly to the parent duration. No displacement is clipped while pretending
that the full time elapsed, and no rejected trial modifies permanent state.
Heun nearly doubles RHS work and neither integrator removes the need for an
internal-step convergence study.

### 10.10 Canonical Geometry, Outlet Priority, and Lifecycle Invariants

Every terminal outlet has an outward unit normal \(\mathbf n_{o,k}\), a point
\(\mathbf x_{o,k}\) on its authoritative section, a tangent, and a finite
aperture. The directed signed distance is

\[
\boxed{
\psi_k(\mathbf x)=(\mathbf x_{o,k}-\mathbf x)\cdot\mathbf n_{o,k}.
}
\]

Thus \(\psi_k>0\) is upstream, \(\psi_k=0\) is on the outlet plane, and
\(\psi_k<0\) is downstream. A segment is an outlet event only when it crosses
the plane in the outward direction and its crossing point lies inside the
section aperture. Labels alone do not terminate a particle.

For each free trial chord, the code independently finds the first solid-wall
event fraction and the first directed outlet fraction. Event ordering is

\[
\boxed{t_o\le t_w\Longrightarrow\text{commit the outlet first}.}
\]

An exact tie at an outlet/side-wall corner therefore belongs to the outlet.
The centre is stored exactly at \(\psi_k=0\), its rotation is advanced only to
that event fraction, the event time is recorded chronologically, and the
permanent ID becomes inactive immediately. The predictive reaction can
redirect a path through an outlet, so the accepted constrained chord is checked
a second time before any domain-membership decision.

Wall distance cannot identify inside versus outside. Every continuing accepted
endpoint must therefore satisfy all three independent conditions:

\[
\boxed{
\mathbf x\in\Omega,
\qquad g_R(\mathbf x)\ge0,
\qquad \mathbf x\in\Omega_R^{\mathrm{in}}.
}
\]

A point outside \(\Omega\) without an authoritative directed outlet crossing
is a path-solver error, not a new valid termination mode. Likewise, an active
particle may never remain downstream with \(\psi_k<0\). The v15 runtime records
outlet events and invariant counters; a valid run has zero active-outside-lumen
and active-outside-accessible-domain violations.

The complete accepted chord is checked analytically against every relevant
finite solid face segment.
If either the endpoint or any interior chord point is non-feasible, the trial is
not projected or clipped. It requests the true-time bisection described in
Section 10.9. Only a small endpoint residual may use the documented exact-
gradient projection; the independently sampled accepted endpoint and complete
chord must then be nonnegative. Exhausting `contact_max_time_refinements` raises an explicit
numerical error.

The immutable finite-face geometry is indexed once when the transport context
is built. Conservative CSR spatial bins feed compiled Numba exact-nearest-face
and swept-disc kernels; the Python implementation remains the reference path.
Ordinary safe free lanes are accepted as one vectorized batch. Only actual
contact, residual projection, outlet, or error lanes enter Python orchestration.

For a contact candidate, the code examines the actual solid faces in
\(\Gamma_w\). Multiple collinear raster faces with one inward normal represent
one wall. Two simultaneous faces with distinct normals represent a corner,
pinch, or opposing wall. Because the hydrodynamic mobility is a single-planar-
wall approximation, this state raises a model-scope error instead of selecting
the nearest face, averaging normals, or assembling an unsupported two-wall
solver. A placeholder direction created for a degenerate distance-field
gradient may keep a bulk RHS finite, but it is never accepted as the normal of
an active rigid-wall constraint; contact without a unique finite normal stops
with an explicit error.

Open faces never belong to \(\Gamma_w\), so an ordinary terminal exit cannot
be misreported as closed-wall or simultaneous-wall contact. Outlet termination
removes the ID from all later RHS and collision evaluations and never reuses its
state slot. New IDs arise only from the precomputed concentration-flux event
schedule; termination itself does not create a replacement. Both schedule
construction and final admission independently require the unregularized inlet
gap to be finite and nonnegative, so an injected state cannot rely on a negative
gap tolerance while the saved invariant reports zero accepted penetrations.

The saved world position remains:

```text
X = origin_x + gx * grid_spacing_um
Y = fixed_y_um
Z = origin_z + gz * grid_spacing_um
```

#### Revised-v15 acceptance invariants

A production trajectory is acceptable only when all of the following hold:

1. Every admitted and time-integrated accepted position has \(g_R\ge0\), and
   `accepted_negative_gap_count == 0`.
2. Each accepted predictive solve has `lambda >= 0` and a small reported
   complementarity residual. A feasible free chord has zero reaction; a
   wall-directed trial uses the closed-form mobility reaction before movement.
3. Every directed outlet event is committed exactly on its plane and wins a
   tie with a wall event. No continuing active particle is outside \(\Omega\),
   outside \(\Omega_R^{\mathrm{in}}\), or downstream of an applicable outlet.
4. A straight frictionless wall preserves nonzero tangential motion. Over
   contact intervals, the cumulative position path and integrated constrained
   velocity path satisfy \(L_{\mathrm{pos}}/L_{\mathrm{vel}}\to1\) as the
   internal step is refined.
5. In free motion, the gap increment converges to the integrated normal
   kinematics, \(\Delta g_R\approx\int\mathbf n\cdot\mathbf V\,dt\). A stable
   positive gap equal to a mobility regularizer must not appear as a geometric
   barrier.
6. Every failed curved chord is retried as two chronological half-time
   intervals. Rejected position, rotation, bond, lifecycle, and accepted-work
   state is never committed.
7. A simultaneous contact with distinct true solid walls fails explicitly as
   outside the single-wall model. It is never downgraded to a warning or
   repaired by arbitrary normal selection.
8. Refining `integration_substeps` and tightening the residual-only
   `contact_geometry_tolerance_um` makes
   along-wall path, target exposure, bond state, and outlet time approach stable
   values.

### 10.11 Runtime Diagnostics and Stability Warning

The console log and trajectory metadata together report the following. The
console prints the most important population, clearance, overlap, timing, and
stability subset; the full aggregate set is retained in metadata:

```text
physical concentration, lambda, and planned injection interval
empty-lumen initial-condition marker and active population from frame 0 onward
unique IDs created and formal admissions/exits
minimum, mean, and maximum active population
inlet wait-event count, current queue, and maximum waiting time
minimum saved wall gap and contact count
minimum accepted internal wall gap and accepted-negative-gap invariant count
contact-constraint evaluations, active reactions, and maximum reaction force
endpoint residual-projection count and maximum projection distance
maximum predictive complementarity residual
physical-time refinement count and maximum depth
contact nonzero-velocity/zero-progress count
contact-only cumulative L_pos, cumulative L_vel, and L_pos/L_vel ratio
maximum free-motion gap-kinematic residual
directed outlet event count
active-outside-lumen and active-outside-Omega_R^in invariant counts
maximum pair physical overlap and collision-layer compression
maximum collision speed and interacting-pair counts
reciprocity mismatch and degenerate near-wall normals
opposite-wall hydrodynamic-validity warnings
particle transport timing
saved-frame interval, internal step, and total internal-step count
maximum internal-step displacement and grid/radius ratios
conservative absolute-displacement/collision-layer ratio
internal-step/collision-relaxation-time ratio
```

For the empty-start public run, the reported minimum active population includes
formal frame 0 and is therefore zero. "Unique bubbles created" counts planned
events whose scheduled times have been reached within the formal horizon; it
can include IDs still waiting upstream and not yet admitted. Maximum inlet wait
is updated only when an ID eventually enters. IDs still waiting at the final
horizon keep NaN admission/wait registry values and do not contribute to that
maximum.

Every saved-frame record re-evaluates the free RHS for output quantities.
Collision and molecular RHS observation counts can therefore include these
readout probes as well as accepted integration stages; they are not counts of
physical intervals. Predictive-contact complementarity, path-length, outlet,
residual-projection, and accepted-gap metrics are accumulated from accepted v15
transactions. The molecular capacity-limited metric is counted only when an
accepted molecular state update is committed. A rejected speculative trial
contributes no accepted position, force, path-length, bond, or elapsed-time
diagnostic; only the administrative refinement request survives.

For accepted contact intervals, the runtime accumulates

\[
L_{\mathrm{pos}}=\sum_n\|\mathbf x_{n+1}-\mathbf x_n\|,
\qquad
L_{\mathrm{vel}}=\sum_n\tau_n\|\mathbf B\mathcal U_n^c\|.
\]

The metadata ratio `contact_position_to_velocity_path_ratio` should approach
one under internal-step refinement. It is computed only where a reaction is
active or the start/end true gap is at roundoff contact; free-flow intervals
do not dilute it. `maximum_free_gap_kinematic_residual_um` separately audits
the discrete free-motion relation
\(\Delta g_R\approx\tau\,\mathbf n\cdot\mathbf V\).

The primary step-size diagnostic compares an accepted v15 constrained motion
proposal with grid spacing and minimum active radius separately. A maximum
ratio at or above `0.2` prints a convergence warning. Contact event splitting
and reaction can redirect the final centre path, so this quantity is a
stage-resolution diagnostic rather than the realized frame velocity.
Absolute proposal displacement divided by the collision layer is retained as a
conservative metadata diagnostic, but it is not the sole warning because
particles can share a large bulk translation without changing their pair
separation. Collision temporal stiffness is reported separately as
`internal_dt_s / collision_relaxation_time_s`, with a warning at `0.2`.

The nominal dynamics subcycle remains configured and does not automatically
adapt to a velocity error estimate. Contact-specific physical-time refinement
is a feasibility mechanism, not a replacement for a global accuracy study.
Production results therefore require explicit
`integration_substeps`, `2*integration_substeps`, and
`4*integration_substeps` convergence studies while keeping output `dt_s` and
`collision_relaxation_time_s` fixed.

### 10.12 Reference Model versus Engineering Extensions

The near-wall mobility and shear-load equations are based on:

```text
references/Optimization of ultrasound contrast agents with computational models to improve.pdf
references/bit_22857_sm_suppdata.pdf
```

The code follows the reduced X-Z translation plus `+Y` rotation formulation
specified in `references/Revised v5.md`. The cited work primarily supplies the
single-sphere, single-planar-wall mobility. Whole-vessel wall/free-space
blending, pair collision relaxation, synchronous integration, deterministic
continuous perfusion, the empty formal start without pre-equilibration, upstream
inlet waiting, the inlet-connected finite-radius domain, directed outlet
lifecycle, opposite-wall validity warnings, and the Revised-v15 predictive
single-wall contact/time-refinement algorithm are explicit engineering
extensions and are identified as such in trajectory metadata.

### 10.13 Revised-v8/v10 Target Selection and Revised-v7 Mean-Field Bond Model

Revised v8 separates candidate generation from biological target assignment.
The program automatically produces **selectable candidate vessel beds**, not an
automatically diagnosed tumour or disease mask. Candidate geometry is therefore
defined by vascular topology and accepted-flow direction; flow fraction,
residence-time scale, and wall shear are descriptive values shown to the user,
not thresholds that silently decide where target molecules exist.

Candidate preparation is a separate run selected by
`--prepare-target-candidates`. It loads and rasterizes the same DCCO vessel
model, solves and validates the steady CFD field, and then stops before cardiac
modulation, particle injection, molecular kinetics, and trajectory rendering.
The accepted vessel graph is oriented by the sign of each saved vessel flow.
Consecutive segments without an inlet, outlet, or branch boundary are grouped
into one maximal basic unit. Numerically zero-flow units remain in the geometry
but do not produce independent perfused candidates.

For every perfused basic unit (e), the catalog contains a local candidate

\[
\mathcal A_e=\{e\},
\]

and, where it adds a larger selectable region, a downstream-subtree candidate

\[
\mathcal S_e=\{e\}\cup\operatorname{Desc}(e).
\]

The complete root subtree and duplicate terminal subtrees are omitted because
they add no useful choice. Ordinary raster cells inherit their existing vessel
ID through the segment-to-unit lookup. Junction-core cells are traced forward
through the accepted velocity field until they reach an unambiguous downstream
unit. A nearest flow-resolved basin fills only those junction cells whose trace
cannot be completed because of numerical stagnation or raster limits; the count
is saved in metadata. Unit ownership is then extended to eligible solid-wall
samples, while open inlet/outlet caps remain excluded.

For a candidate \(\mathcal C\), the catalog reports

\[
Q_{\mathcal C}=|Q_e|,
\qquad
\phi_{\mathcal C}=\frac{|Q_e|}{|Q_{\mathrm{root}}|},
\qquad
V_{\mathcal C}=\sum_{f\in\mathcal C}\pi r_f^2L_f,
\qquad
T_{\mathcal C}=\frac{V_{\mathcal C}}{Q_{\mathcal C}},
\]

plus a wall-area-weighted mean wall-shear proxy. The cylindrical endothelial
surface estimate is

\[
A_{w,e}=2\pi r_eL_e,
\qquad
\alpha_j=\frac{\sum_{e\in\mathcal C_j}A_{w,e}}
{\sum_{e\in E}A_{w,e}}.
\]

The denominator contains every network segment, including unperfused geometry.
Schema v3 additionally distributes each segment area \(2\pi r_eL_e\) over its
associated closed raster-wall sites. If \(w_i\) is the resulting quadrature
weight of wall site \(i\), area fractions are evaluated from \(\sum_i w_i\),
not from the number of wall pixels. Segment area that cannot be associated with
its own wall sites is first redistributed within the same natural vessel unit;
any remaining unmapped theoretical area is reported explicitly.

For multiple root vessels, the network inlet reference is the sum of all root
flows. The expected number of injected bubbles entering unit or candidate
\(j\) is

\[
N_j^{\mathrm{exp}}=\lambda_{\mathrm{in}}T_{\mathrm{obs}}
\frac{Q_j}{\sum_r Q_{\mathrm{root},r}},
\qquad
T_{\mathrm{obs}}=n_{\mathrm{steps}}\,\Delta t_{\mathrm{output}}.
\]

Here \(\lambda_{\mathrm{in}}\) is rebuilt with the same finite-size inlet
cross-section model used by formal continuous perfusion. A wall site is
analytically accessible only when its natural vessel unit is perfused and has
\(N_j^{\mathrm{exp}}\geq1\). Single-vessel candidates remain available for
manual selection.

Revised v10 uses a downstream subtree only as a **coarse location anchor**. It
is never interpreted as the final all-positive target. Given the requested
influence-region network wall-area fraction \(\alpha_G\), eligible anchors are
ranked lexicographically, without a weighted empirical score:

\[
\operatorname*{arg\,min}_j
\left(
\left|\ln\frac{\alpha_j}{\alpha_G}\right|,
-d_j,
R_{g,j},
\operatorname{ID}_j
\right).
\]

The wall-area radius of gyration \(R_g\) is calculated from segment centres and
the same \(2\pi rL\) weights. Size match is decisive first, followed only on
exact ties by deeper topology, smaller spatial spread, and stable ID. Flow is
used only for the analytical accessibility guard; speed, residence time, and
WSS do not determine target-positive probability.

The selected anchor's wall-area centroid gives the physical X-Z centre
\(\mathbf c_G\) of a circular tissue-space influence region. Its radius is the
available wall-site distance whose cumulative physical weight most closely
matches \(\alpha_G A_{w,\mathrm{network}}\). Therefore the circle can intersect
several vessels that are close in tissue space even when they are unrelated in
the topology. Only accessible closed wall inside the circle is eligible:

\[
\Omega_G=\{\mathbf x:\|\mathbf x-\mathbf c_G\|\le R_G\},
\qquad
\Gamma_G=\Gamma_{w,\mathrm{accessible}}\cap\Omega_G.
\]

On physical wall coordinates in \(\Gamma_G\), the implementation evaluates a
seeded Random Fourier Feature approximation to a zero-mean Gaussian field with
squared-exponential covariance

\[
C(\mathbf x,\mathbf x')=
\exp\!\left(-\frac{\|\mathbf x-\mathbf x'\|^2}{2\ell_T^2}\right).
\]

The realization is thresholded so that the requested fraction \(\phi_T\) is
matched by **physical wall area**:

\[
\frac{\sum_{i\in\Gamma_G}w_i I[Z_i\ge z_T]}
     {\sum_{i\in\Gamma_G}w_i}
\approx\phi_T,
\qquad
\chi_T(i)=I[i\in\Gamma_G]I[Z_i\ge z_T].
\]

The nearest attainable weighted fraction is used because one indivisible wall
site can straddle the requested area threshold. The correlation length
\(\ell_T\) is expressed in micrometers. The PCG64 seed and Random Fourier
wavevectors/phases are saved, so identical geometry, grid, parameters, seed,
and mode count reproduce the same realization. Multiple seeds still need to be
run explicitly for biological-layout sensitivity analysis.

These rules generate a repeatable **spatially heterogeneous synthetic molecular
target**, not a tumour diagnosis. Automatic mode writes both the canonical
Boolean target NPZ and a schema-v2 audit JSON. The Trame interface has explicit
manual and automatic workflows: automatic preview highlights the physical
influence region, accessible influence walls, and final positive patches;
checking the candidate tree switches to the independent manual workflow.

In manual mode, the descriptive values help the user understand a choice but do
not change its membership. The Trame selector uses server-side VTK rendering so
large planar grids do not depend on client-side mesh serialization. It displays
the hierarchy over the accepted speed field with sparse flow arrows, removes
child selections already contained by a selected parent, and exports the
Boolean union

\[
M_{\mathrm{target}}(x,z)=\bigvee_{c\in\mathcal I_{\mathrm{sel}}}M_c(x,z).
\]

The formal simulation accepts only `mask_npz` (or disabled `none`). Its Boolean
ROI must be accompanied by strictly increasing physical X and Z coordinate axes
in micrometers. The selector writes the canonical keys `x_um`, `z_um`, and
`target_mask`; externally justified masks may use configured aliases but must
obey the same physical-coordinate contract. When the saved axes exactly match
the current CFD grid, target cells are mapped directly. A mask on another
physical grid is sampled at the registered grid overlay by nearest neighbour.
The former `physical_polygon` mode
and `polygon_vertices_um` input have been removed and are rejected by the
configuration loader.

The optional Revised-v7 extension then separates **where target molecules
exist** from **how bonds form and transmit load**.

The input ROI does not directly become a force-producing image region. It is
intersected with the canonical finite face segments in \(\Gamma_w\). Open CFD
faces have already been removed from \(\Gamma_w\), so no separate opening-mask
dilation is applied and legitimate side walls next to an outlet remain
eligible. Consequently,

\[
\mathcal T \subseteq \mathcal W_{\mathrm{solid}},
\qquad
\mathcal T \cap \mathcal O_{\mathrm{open}} = \varnothing ,
\]

where \(\mathcal T\) is the target-positive subset of the authoritative solid
faces. The saved cell mask remains the selector/visualizer overlay, while a
separate Boolean flag is constructed once for every solid face so runtime
queries never resolve a half-grid nearest-cell tie differently on opposite
walls. Saved eligible wall coordinates are the physical centres of those face
segments; their axis, full length, inward normal, and target-positive flag are
saved alongside them. The density is fixed at the configured
\(\rho_T\) in `molecule/m^2` on selected sites and zero elsewhere; targets are
not depleted as bubbles bind.

For a bubble of radius \(R\), molecular capture uses the same true gap
\(g_R=d_w-R\) as particle contact. A materially negative gap is rejected; only
floating-point roundoff at zero may be clamped. With capture distance \(d_c\),
the local overlap depth and reaction-disk radius are

\[
\delta=\max(d_c-\max(g_R,0),0),
\qquad
a=\sqrt{\max(2R\delta-\delta^2,0)}.
\]

The runtime requires \(d_c\le 2R_i\) for every evaluated bubble. In a
polydisperse run, the configured capture distance must therefore not exceed the
smallest scheduled/evaluated bubble diameter. Before target evaluation, the
particle centre is projected directly onto the nearest finite segment of
\(\Gamma_w\); no wall point is first reconstructed from a bilinear distance and
gradient. At an exact equal-distance corner, every tied canonical face is
evaluated independently and the largest target-positive reaction area is used.
This deterministic rule keeps the overlap fraction at or below one and never
averages distinct wall normals; a true simultaneous rigid contact remains a
hard model error in the contact solver. Let \(\chi_T\) be the per-face target
indicator along each candidate canonical wall tangent. The
target-positive portion of the reaction disk is found from the exact union of
target-positive closed-wall intervals and the resulting circular segments are
integrated analytically:

\[
A_{\mathrm{rxn}}
=2\int_{-a}^{a}\chi_T(s+\eta)\sqrt{a^2-\eta^2}\,d\eta .
\]

For a completely target-positive local plane, this gives
\(A_{\mathrm{rxn}}=\pi a^2\). The code computes this area in `um^2` for output
and converts it to `m^2` before combining it with molecular surface densities.
"Analytic" here means that there is no residual fixed-node quadrature error for
the given interval representation. Eligible tangent support is the exact union
of collinear canonical solid face segments within the reaction disk; there is
no `0.8*grid_spacing_um` proximity band and no two-cell opening dilation. The
source Boolean ROI remains grid-sampled, and \(\Gamma_w\) itself comes from the
raster CFD boundary, so reaction-area convergence must still be checked against
`domain.grid_spacing_um`.
The expected local ligand and target populations before subtracting existing
bonds are

\[
N_L=\rho_L A_{\mathrm{rxn}},
\qquad
N_T=\rho_T A_{\mathrm{rxn}},
\]

where the area in these products is expressed in square metres.

The bond state is a deterministic two-moment mean field, not a list of discrete
bonds:

- \(n\ge 0\) is the non-integer expected bond count; and
- \(m\) is the signed total tangential bond extension in micrometers.

When \(n>0\), the mean tangential extension is \(q=m/n\). When \(n=0\), the
canonical state is \(n=m=q=0\). With the effective two-surface association
coefficient \(k_{\mathrm{on}}^{2D}\) in `m^2/(molecule s)`, bond formation is

\[
r_+
=\frac{k_{\mathrm{on}}^{2D}}{A_{\mathrm{rxn}}}
  [N_L-n]_+[N_T-n]_+,
\]

and is zero when the reaction area or either available population is zero. A
paper value reported only in `s^-1` is not dimensionally interchangeable with
this coefficient.

The current bond length, tensile extension, single-bond tension, and Bell
dissociation rate are

\[
L=\sqrt{\bar g^2+q^2},
\qquad
e=[L-\ell_0]_+,
\qquad
f=k_b e,
\]

\[
k_{\mathrm{off}}(f)
=k_{\mathrm{off}}^0
 \exp\!\left(\frac{f x_b}{k_B T}\right).
\]

Here \(k_b\) is in `pN/um`, \(x_b\) is configured in nanometres, and the code
applies the `pN nm` to joule conversion before forming the dimensionless Bell
exponent. The configured exponent cap is a numerical overflow guard and its use
is reported diagnostically.

Let \(V_t\) be particle translation along the local tangent and
\(\Omega_y\) its angular velocity about `+Y`. The surface slip used to stretch
the bond population is

\[
v_{\mathrm{slip}}=V_t+R\Omega_y.
\]

The two state equations are

\[
\frac{dn}{dt}=r_+-k_{\mathrm{off}}n,
\qquad
\frac{dm}{dt}=n v_{\mathrm{slip}}-k_{\mathrm{off}}m.
\]

The resulting local tangential force, inward-normal component, and torque are

\[
F_t=-n f\frac{q}{L},
\qquad
F_n=-n f\frac{\bar g}{L},
\qquad
T_y=R F_t.
\]

Because the saved wall normal points from wall into lumen, negative \(F_n\)
pulls the particle toward the wall. The global X-Z molecular force is
\(\mathbf F_b=F_t\mathbf t+F_n\mathbf n\).

The molecular force and torque are added to the collision load and passed
through the same position-, radius-, viscosity-, and gap-dependent mobility
operator. Thus bonds alter both translation and rotation, but do not modify or
re-solve the one-way-coupled CFD field. The bond right-hand side is first
evaluated with the preliminary slip, the mobility update is applied, and the
pure bond kinetics are then re-evaluated with the resulting particle slip before
the synchronous state update. Exponential Euler or exponential Heun integrates
the linear dissociation sink while preserving nonnegative \(n\). The accepted
contact-constrained tangential displacement and accepted rotation, not a free
or rejected trial, provide the extension source for \(m\). The predictor is
made feasible and its stage-2 RHS is constrained before the Heun correction.
Molecular state is committed only after the final geometry and physical-time
interval have been accepted.

For this accepted-motion extension update, Euler uses the start-state wall
tangent. Heun normalizes the sum of start and predictor normals to obtain one
representative tangent, and uses the mean of old and predicted expected bond
counts. Centre translation is the accepted predictive-step displacement;
rotation is the corresponding accepted angular displacement. A directed-outlet
split or physical-time refinement integrates bond kinetics only over accepted
chronological pieces. Position-root iterations, chord checks, and rejected
full-step trials add no molecular capture time and do not update \(n\) or
\(m\). No third RHS is evaluated at the final accepted position.

The formation-capacity limiter uses
`max(old expected bond count, instantaneous capacity)`: it can limit newly
formed \(n\), but it does not delete old bonds or proportionally rescale the
already integrated \(m\) when reaction area shrinks. Formation then stops and
dissociation continues. Every newly admitted permanent ID starts with
\(n=m=0\), and terminated state slots are never reused.

Only the two moments \(n,m\) are persistent. There is no discrete receptor
identity or fixed wall-anchor coordinate. As a bubble moves, existing mean
bonds are reinterpreted using the current nearest-wall normal, tangent, and
gap. This is a material limitation near corners, bifurcations, or a nearest-wall
switch and must not be described as explicit tracking of individual anchored
bonds.

### 10.14 Empty-Lumen No-Bond Contact Pilot and `Da_on` Scenarios

When experimental two-dimensional association data are unavailable, the
implemented alternative is a clearly labelled sensitivity workflow rather than
an invented calibrated rate. The pilot uses the normal empty-lumen continuous
perfusion run with `molecular_binding.enabled=false`, so the target geometry is
present but exerts no molecular force. Full particle diagnostics are required
because the pilot reconstructs tangential surface slip from translation and
angular velocity.

The pilot uses the saved same-permanent-ID realized centre velocity and the
saved instantaneous angular velocity,
\(v_{\mathrm{slip}}=V_{t,\mathrm{realized}}+R\Omega_y\), at each accepted
saved position. The centre term therefore includes the effect of wall geometry
acceptance instead of using a non-realized mobility request. A degenerate saved
wall normal is excluded from pilot contact by forcing its reaction radius to
zero. In the formal binding RHS, by contrast, a degenerate distance-field
normal uses the counted deterministic `+Z` fallback basis; the pilot rule must
not be assumed to apply to formal binding.

For each configured ratio \(d_c/\ell_0\), a saved active observation is counted
as target contact only when its target-positive reaction area is greater than
zero. Each positive saved observation contributes one formal output interval
`particles.dt_s` by an explicit rectangle rule. Contiguous positive samples form
contact events and per-bubble cumulative exposure. Reaction area and slip are
reported as contact-time-weighted means. An event that is still positive in the
final saved frame is explicitly marked right-censored, and the report gives its
count, duration, and fraction of total observed contact time. These quantities
describe transport exposure only; none is used to define a molecular rate.

`binding_scenario_sweep.da_on_reference_time_s` must be fixed before transport.
For each requested dimensionless level, target-density level, ligand-density
level, and capture/rest-length ratio, the scenario builder uses

\[
k_{\mathrm{on}}^{2D}
=\frac{Da_{\mathrm{on}}}
       {\rho_T t_{\mathrm{ref}}},
\qquad
t_{\mathrm{ref}}=\texttt{da\_on\_reference\_time\_s}.
\]

Target-density sweep values entered in `molecule/um^2` are converted using

\[
1\ \mathrm{molecule/um^2}=10^{12}\ \mathrm{molecule/m^2}.
\]

The configured ligand-density levels are retained as a separate sensitivity
axis. Scenario construction reads only these predeclared inputs and never reads
the pilot's observed contact time; zero observed contact therefore remains a
valid exposure result and does not alter the scenario table. The fixed
\(t_{\mathrm{ref}}\) is a declared dimensionless scaling choice, not a fitted
residence time or evidence of calibrated biology.

Numerical wall lock remains a secondary exposure-quality audit. For current
v15 trajectories it compares the last accepted internal constrained velocity
with same-ID realized saved-frame motion when the wall reaction is active.
Historical v13 files retain their read-only legacy predicate. The primary v14
transport audit is the internal contact
\(L_{\mathrm{pos}}/L_{\mathrm{vel}}\) diagnostic, which is not limited by
saved-frame cadence. The resulting YAML files are definitions and exposure
reports for hypothesis generation; they neither fit a desired adhesion outcome
nor automatically execute bound-particle simulations. Contact-duration results
still require output-cadence convergence, independently of internal-step
convergence.

### 10.15 Revised-v7/v8 Reproduction Checks

For every molecular run, archive and report all of the following:

1. The candidate catalog, selected candidate IDs, final target-mask NPZ, its
   physical coordinate axes, and the external reasoning used to assign the
   selected candidate vessel beds as biologically target-positive.
2. The configured target and ligand densities, capture distance, rest length,
   effective two-dimensional association coefficient, zero-force off-rate,
   stiffness, reactive compliance, temperature, and the provenance and units of
   every value.
3. The number of saved target wall sites, confirmation that inlet/outlet caps
   are excluded, and visual or numerical agreement between
   `molecular_target_field.npz` and the intended wall region.
4. The minimum expected ligand and target counts in contacted reaction areas.
   Values below `mean_field_warning_count` do not disable binding, but they
   weaken the deterministic mean-field approximation and should trigger a
   sensitivity analysis or a future stochastic model.
5. Fixed-output-step and internal-substep convergence checks for position,
   expected bond count, total extension, force, torque, and contact duration.
6. Whether the Bell exponent cap or formation-capacity limiter was reached and
   whether capture distance or reaction-disk size violated the reported local
   planar-wall validity diagnostics.
7. For a no-experiment study, declare `da_on_reference_time_s`, `Da_on`, and
   density axes before transport. Use the no-bond pilot only to quantify target
   exposure; do not replace the declared reference time with its observed
   contact duration. Keep every scenario labelled as sensitivity input rather
   than a calibrated biological constant.

## 11. Output NPZ Contents

### 11.1 `velocity_and_wall_shear_field.npz`

| Field                                                                                         |            Shape | Meaning                                                                                                                                                                                              |
| --------------------------------------------------------------------------------------------- | ---------------: | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `field_schema_version`                                                                      |           scalar | `v15_exact_solid_face_particle_boundary` for a current accepted field containing all authoritative particle boundary arrays.                                                                       |
| `origin_um`                                                                                 |         `(3,)` | World coordinate of the grid origin.                                                                                                                                                                 |
| `spacing_um`                                                                                |         `(1,)` | Grid spacing in micrometers.                                                                                                                                                                         |
| `fixed_y_um`                                                                                |         `(1,)` | Constant Y coordinate used by the 2D X-Z model.                                                                                                                                                      |
| `shape`                                                                                     |         `(2,)` | Grid shape`(nx, nz)`.                                                                                                                                                                              |
| `x_coordinates_um`                                                                          |        `(nx,)` | World X coordinate of grid cell centers.                                                                                                                                                             |
| `z_coordinates_um`                                                                          |        `(nz,)` | World Z coordinate of grid cell centers.                                                                                                                                                             |
| `lumen_mask`                                                                                |     `(nx, nz)` | True for rasterized vessel-lumen cells.                                                                                                                                                              |
| `lumen_fraction`                                                                            |     `(nx, nz)` | Exact Shapely area fraction \(\operatorname{Area}(C_{ij}\cap\Omega)/\operatorname{Area}(C_{ij})\).                                                                                                 |
| `wall_mask`                                                                                 |     `(nx, nz)` | Grid cells near the lumen boundary.                                                                                                                                                                  |
| `junction_core_mask`                                                                        |     `(nx, nz)` | Junction-core cells generated by continuous geometry, when available.                                                                                                                                |
| `vessel_id`                                                                                 |     `(nx, nz)` | Zero-based nearest vessel ID for each lumen cell,`-1` outside. Trajectory `record_vessel_id` intentionally uses a separate 1-based convention.                                                   |
| `radius_um`                                                                                 |     `(nx, nz)` | Local DCCO vessel radius assigned during rasterization.                                                                                                                                              |
| `flow_rate_um3_s`                                                                           |     `(nx, nz)` | Local vessel flow rate assigned during rasterization.                                                                                                                                                |
| `q2d_flow_um2_s`                                                                            |     `(nx, nz)` | Local 2D flow target`flow_rate / effective_thickness`.                                                                                                                                             |
| `viscosity_mpas`                                                                            |     `(nx, nz)` | Local apparent viscosity assigned during rasterization.                                                                                                                                              |
| `distance_to_wall_um`                                                                       |     `(nx, nz)` | Exact distance from each accepted grid centre to the continuous solid wall in micrometers; zero outside the Boolean lumen mask.                                                                     |
| `wall_normal_xz`                                                                            |  `(nx, nz, 2)` | Inward normal sampled from the continuous solid wall at each accepted grid centre.                                                                                                                  |
| `velocity_xz_um_s`                                                                          |  `(nx, nz, 2)` | Accepted cell-centered X-Z velocity reconstructed from final face fluxes.                                                                                                                            |
| `speed_um_s`                                                                                |     `(nx, nz)` | Velocity magnitude.                                                                                                                                                                                  |
| `initial_velocity_xz_um_s`                                                                  |  `(nx, nz, 2)` | Initialized pre-iteration X-Z velocity, saved when available.                                                                                                                                        |
| `initial_speed_um_s`                                                                        |     `(nx, nz)` | Magnitude of the initialized velocity, saved when available.                                                                                                                                         |
| `wall_shear_stress_pa`                                                                      |     `(nx, nz)` | In-plane WSS proxy.                                                                                                                                                                                  |
| `divergence_s_inv`                                                                          |     `(nx, nz)` | Final finite-volume divergence from saved face fluxes.                                                                                                                                               |
| `wall_penetration_um_s`                                                                     |     `(nx, nz)` | Solid-wall normal velocity diagnostic; open boundaries excluded.                                                                                                                                     |
| `pressure`                                                                                  |     `(nx, nz)` | Pressure-like projection potential in arbitrary units, not a calibrated physical pressure in Pa.                                                                                                     |
| `inlet_label`, `outlet_label`                                                             |     `(nx, nz)` | Labeled cells adjacent to selected open boundary faces.                                                                                                                                              |
| `boundary_velocity_xz_um_s`                                                                 |  `(nx, nz, 2)` | Open-boundary face flux represented as a cell-centered helper field.                                                                                                                                 |
| `boundary_normal_xz`                                                                        |  `(nx, nz, 2)` | Accumulated open-boundary normals.                                                                                                                                                                   |
| `boundary_weight`                                                                           |     `(nx, nz)` | Boundary profile weights.                                                                                                                                                                            |
| `boundary_edge_length_um`                                                                   |     `(nx, nz)` | Accumulated selected open-face lengths.                                                                                                                                                              |
| `open_boundary_flux_um2_s`                                                                  |     `(nx, nz)` | Signed open-boundary flux accumulated per boundary cell.                                                                                                                                             |
| `face_flux_x_um2_s`                                                                         | `(nx + 1, nz)` | Final x-face fluxes.                                                                                                                                                                                 |
| `face_flux_z_um2_s`                                                                         | `(nx, nz + 1)` | Final z-face fluxes.                                                                                                                                                                                 |
| `inlet_target_by_label_um2_s`, `outlet_target_by_label_um2_s`                             |         variable | Per-label target fluxes.                                                                                                                                                                             |
| `inlet_actual_by_label_um2_s`, `outlet_actual_by_label_um2_s`                             |         variable | Per-label actual final fluxes integrated from face fluxes.                                                                                                                                           |
| `open_face_cell_ij`, `open_face_index_ij`                                                 |        `(B,2)` | Lumen-side cells and finite-volume face indices for every authoritative open face.                                                                                                                   |
| `open_face_axis`, `open_face_label`, `open_face_kind`                                   |         `(B,)` | Face axis, boundary label, and kind (`-1` inlet, `+1` outlet). These—not raster label bands—define which boundary faces are removed from \(\Gamma_w\).                                         |
| `open_face_normal_xz`, `open_face_center_xz_um`                                           |        `(B,2)` | Outward unit normal and physical face centre for every authoritative open face.                                                                                                                      |
| `open_face_length_um`                                                                       |         `(B,)` | Physical measure of each selected open face.                                                                                                                                                         |
| `open_section_point_xz_um`, `open_section_outward_normal_xz`, `open_section_tangent_xz` |        `(S,2)` | Directed inlet/outlet section geometry used by particle lifecycle.                                                                                                                                   |
| `open_section_half_width_um`, `open_section_label`, `open_section_kind`                 |         `(S,)` | Section aperture, label, and kind. Outlet rows define the\(\psi_k\) tests.                                                                                                                           |
| `solver_<key>`                                                                              |  usually`(1,)` | Solver metadata wrapped with`np.asarray([value])`; list-valued metadata may therefore have additional dimensions.                                                                                  |

`domain_metadata.yaml` additionally records
`particle_boundary_geometry_schema=v16_continuous_swept_vessel_boundary`,
continuous wall-element counts, open-section counts, curve tessellation
settings, and the continuous geometry hash. There is no raster-wall fallback
or legacy geometry mode.

### 11.2 `microbubble_field_trajectories.npz`

The current writer uses flat append-only observation records plus a permanent-ID
lifecycle registry. Define:

```text
F   = n_steps + 1               stored frames, including formal time zero
N_f = active bubbles in frame f
R   = sum_f N_f                 flat active-observation records
U   = scheduled permanent IDs through the formal horizon
```

`frame_offsets[f]:frame_offsets[f+1]` is the flat slice for frame `f`. Frame
population is not prescribed: it emerges from the concentration-derived inlet
rate, residence time, outlet exits, and upstream waiting. Every stored record is
an active bubble in that frame. An event that cannot enter without physical
overlap remains in a deterministic ordered-retry queue; a later unblocked event
may be admitted while an earlier event remains blocked, so admission is not
strict FIFO. No overlapping fallback is inserted into the lumen.

#### Required arrays

| Field                                               |        Shape | Meaning                                                                                                                                                                                                                                                                                                           |
| --------------------------------------------------- | -----------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `frame_offsets`                                   | `(F + 1,)` | Monotone int64 offsets into all`record_*` arrays.                                                                                                                                                                                                                                                               |
| `record_bubble_id`                                |     `(R,)` | Permanent int64 bubble ID assigned by deterministic injection order; an ID and its state slot are never reused.                                                                                                                                                                                                   |
| `record_positions_um`                             |   `(R, 3)` | Bubble-center XYZ world position; Y is fixed by the planar model.                                                                                                                                                                                                                                                 |
| `record_velocities_um_s`                          |   `(R, 3)` | Last accepted internal v15 constrained generalized velocity for the permanent ID; a newly admitted ID with no accepted step uses the current saved-state RHS fallback. It includes carrier/wall-shear, collision, optional bond, and unilateral reaction, but is not the realized frame-to-frame centre velocity. |
| `record_realized_velocities_um_s`                 |   `(R, 3)` | Same-permanent-ID frame-to-frame centre velocity derived from geometry-accepted positions. It is the centre-translation input used by the no-bond contact-pilot slip calculation.                                                                                                                                 |
| `record_wall_shear_stress_pa`                     |     `(R,)` | Bilinear sample of the nonnegative in-plane WSS proxy, multiplied by the local cardiac factor when enabled.                                                                                                                                                                                                       |
| `record_vessel_id`                                |     `(R,)` | Nearest raster vessel ID stored as 1-based for the active observation.                                                                                                                                                                                                                                            |
| `record_active`                                   |     `(R,)` | True for every stored continuous-perfusion record; retained for schema compatibility.                                                                                                                                                                                                                             |
| `record_diameter_um`                              |     `(R,)` | Immutable diameter belonging to`record_bubble_id`.                                                                                                                                                                                                                                                              |
| `record_wall_gap_um`                              |     `(R,)` | Canonical bilinear sample of the oriented signed wall-distance field minus particle radius,\(g_R=d_w-R\). Its zero-distance contour is authoritative \(\Gamma_w\), and every current accepted active record must have nonnegative gap.                                                                            |
| `record_wall_contact`                             |     `(R,)` | Active and`wall_gap <= wall_contact_threshold_um`. It is a diagnostic, not a hard safety margin.                                                                                                                                                                                                                |
| `record_wall_normal_xz`                           |   `(R, 2)` | Inward normal from nearest solid wall into lumen; opposite to field-NPZ raster normal convention.                                                                                                                                                                                                                 |
| `record_contact_constraint_active`                |     `(R,)` | True when the last accepted internal predictive solve used a strictly positive reaction. It is distinct from the broader output-only`record_wall_contact` proximity flag.                                                                                                                                       |
| `record_contact_reaction_force_pn`                |     `(R,)` | Last accepted internal position reaction magnitude`lambda` in pN. It is nonnegative and zero for a feasible free chord.                                                                                                                                                                                         |
| `record_contact_free_normal_velocity_um_s`        |     `(R,)` | `a.T @ U0` from the last accepted internal interval before adding the rigid-wall reaction.                                                                                                                                                                                                                      |
| `record_contact_constrained_normal_velocity_um_s` |     `(R,)` | `a.T @ U^c` for the last accepted internal v15 interval. It is a diagnostic of that accepted velocity, while endpoint complementarity is audited separately in metadata.                                                                                                                                        |
| `registry_bubble_id`                              |     `(U,)` | Permanent IDs`0..U-1`.                                                                                                                                                                                                                                                                                          |
| `registry_diameter_um`                            |     `(U,)` | One immutable diameter per permanent ID.                                                                                                                                                                                                                                                                          |
| `birth_frame`, `death_frame`                    |     `(U,)` | Lifecycle frame indices relative to empty formal frame 0. Registry arrays below provide the finer scheduled/admission/path-event times.                                                                                                                                                                           |
| `termination_reason`                              |     `(U,)` | Current uint8 codes are`0=none` and `1=authoritative directed outlet`. A non-outlet domain departure is a v15 numerical error rather than a committed lifecycle reason; historical readers may still encounter archived code `3`.                                                                           |
| `active_count_per_frame`                          |     `(F,)` | Active observation count for every stored frame.                                                                                                                                                                                                                                                                  |
| `injected_count_per_frame`                        |     `(F,)` | Events admitted during the preceding formal output interval; frame 0 is zero because the formal simulation starts empty.                                                                                                                                                                                          |
| `terminated_count_per_frame`                      |     `(F,)` | Permanent IDs terminated at the transition into that frame.                                                                                                                                                                                                                                                       |
| `registry_scheduled_injection_time_s`             |     `(U,)` | Deterministic scheduled event time relative to the empty-lumen formal time zero; values are nonnegative.                                                                                                                                                                                                          |
| `registry_admission_time_s`                       |     `(U,)` | Lumen-admission time at a scheduled event or accepted segment boundary after waiting; NaN if not admitted by the formal horizon.                                                                                                                                                                                  |
| `registry_exit_time_s`                            |     `(U,)` | Chronological authoritative outlet-plane event time; NaN while alive or not admitted. Historical files may contain a time paired with archived non-outlet reason`3`.                                                                                                                                            |
| `registry_inlet_wait_time_s`                      |     `(U,)` | Nonnegative admission minus scheduled time; NaN until admission, so final still-waiting IDs have no finite wait duration here.                                                                                                                                                                                    |
| `metadata_keys`, `metadata_values`              |     variable | Unicode key/value arrays. Metadata values are stringified during NPZ writing.                                                                                                                                                                                                                                     |

#### Cardiac arrays

These optional keys exist whenever `cardiac_pulsatility.enabled=true` and are
independent of `store_full_diagnostics`:

| Field                           |      Shape | Meaning                                                                                   |
| ------------------------------- | ---------: | ----------------------------------------------------------------------------------------- |
| `record_cardiac_multiplier`   |   `(R,)` | Local retarded-phase multiplier actually used for each stored active observation.         |
| `cardiac_waveform_time_s`     |   `(C,)` | Time coordinates for one compact cardiac cycle.                                           |
| `cardiac_waveform_multiplier` |   `(C,)` | Positive periodic, normally cycle-mean-one multiplier samples.                            |
| `cardiac_path_distance_um`    | grid shape | Root-to-grid vascular path distance used for propagation.                                 |
| `cardiac_delay_s`             | grid shape | Physical propagation delay`cardiac_path_distance_um / pulse_propagation_velocity_um_s`. |

#### Optional full mobility diagnostics

These keys exist only for mobility transport with
`particle_dynamics.store_full_diagnostics=true`:

| Field                                                  |      Shape | Meaning                                                                                                                                                                                               |
| ------------------------------------------------------ | ---------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `record_fluid_velocities_um_s`                       | `(R, 3)` | Local carrier-fluid velocity after cardiac modulation, separate from the particle mobility RHS.                                                                                                       |
| `record_angular_velocity_rad_s`                      |   `(R,)` | Instantaneous constrained particle angular RHS about`+Y`, including near-wall mobility, collision, optional molecular force/torque, and any translation-rotation effect of the unilateral reaction. |
| `record_rotation_angle_rad`                          |   `(R,)` | Integrated, unwrapped rotation angle for each active ragged observation.                                                                                                                              |
| `record_collision_force_xz_pn`                       | `(R, 2)` | Total synchronous pair-repulsion force in the X-Z plane, in pN.                                                                                                                                       |
| `record_collision_neighbor_count`                    |   `(R,)` | Number of interacting pair neighbors in this RHS state.                                                                                                                                               |
| `record_gap_ratio`                                   |   `(R,)` | Signed raw wall gap divided by bubble radius.                                                                                                                                                         |
| `record_near_wall_weight`                            |   `(R,)` | Wall/free-space blend weight:`1` is near-wall, `0` is bulk.                                                                                                                                       |
| `record_opposite_wall_hydrodynamic_validity_warning` |   `(R,)` | Diagnostic that the nearest-single-planar-wall hydrodynamic approximation may be inadequate; this is distinct from true simultaneous wall contact, which stops the run.                               |

#### Molecular bond arrays

These keys exist whenever `molecular_binding.enabled=true` and are independent
of `particle_dynamics.store_full_diagnostics`:

| Field                                                 |      Shape | Meaning                                                                                                                                          |
| ----------------------------------------------------- | ---------: | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `record_bond_count_expected`                        |   `(R,)` | Deterministic expected bond count\(n\); it is continuous and must not be rounded to a binary or integer state.                                   |
| `record_bond_total_tangential_extension_um`         |   `(R,)` | Signed total tangential extension\(m\) in micrometers.                                                                                           |
| `record_bond_mean_tangential_extension_um`          |   `(R,)` | Mean extension\(q=m/n\), canonicalized to zero when no expected bonds exist.                                                                     |
| `record_bond_force_xz_pn`                           | `(R, 2)` | Total molecular force in the global X-Z plane in pN.                                                                                             |
| `record_bond_force_tangent_pn`                      |   `(R,)` | Local wall-tangential molecular force in pN.                                                                                                     |
| `record_bond_force_normal_pn`                       |   `(R,)` | Molecular-force component along the saved inward normal in pN; a negative value pulls toward the wall.                                           |
| `record_bond_torque_pn_um`                          |   `(R,)` | Molecular torque about`+Y` in `pN um`.                                                                                                       |
| `record_single_bond_tension_pn`                     |   `(R,)` | Current Hookean tension of one representative mean bond in pN.                                                                                   |
| `record_bond_formation_rate_bonds_s`                |   `(R,)` | Realized bond-formation rate\(r_+\) in expected bonds per second; this is not the configured `association_rate_m2_per_molecule_s` coefficient. |
| `record_bond_dissociation_rate_s_inv`               |   `(R,)` | Bell-modified first-order dissociation rate in`s^-1`.                                                                                          |
| `record_target_reaction_area_um2`                   |   `(R,)` | Target-positive part of the local reaction disk in`um^2`.                                                                                      |
| `record_available_ligand_count`                     |   `(R,)` | Expected ligand population still available to bond formation.                                                                                    |
| `record_available_target_count`                     |   `(R,)` | Expected target population still available to bond formation.                                                                                    |
| `record_target_overlap_fraction`                    |   `(R,)` | Fraction of the geometrically available reaction disk that is target-positive.                                                                   |
| `registry_final_bond_count_expected`                |   `(U,)` | Final expected bond count stored for every scheduled permanent ID.                                                                               |
| `registry_final_bond_total_tangential_extension_um` |   `(U,)` | Final signed total extension stored for every scheduled permanent ID.                                                                            |

The current schema identifier is:

```text
mobility continuous perfusion without molecular binding:
    continuous_perfusion_empty_start_records_v9
mobility continuous perfusion with molecular binding:
    continuous_perfusion_empty_start_molecular_records_v10
```

The continuous-perfusion metadata records the effective time plan, inlet
concentration and flux conversion, population and waiting statistics, selected
acceleration backend, numeric-kernel family, Numba version/thread capacity,
total particle-transport time, internal-step throughput, selected collision
search paths, wall/contact and collision diagnostics, proposal-displacement
resolution ratios, velocity semantics, engineering-extension labels, and
`state_storage_slots_reused=false`. The selected backend records
`particle_numeric_kernel_family=numba_batched_component_kernels_v18` when
Numba is active, or `python_reference` for the reference path. Revised-v15
geometry, contact, kinematic, and lifecycle auditing includes the following
fixed metadata keys:

```text
wall_contact_integrator = revised_v15_predictive_mobility_unilateral_single_wall
maximum_simultaneous_wall_constraints = 1
```

```text
wall_contact_integrator
maximum_simultaneous_wall_constraints
contact_geometry_tolerance_um
contact_constraint_evaluations
active_contact_constraint_evaluations
maximum_contact_reaction_force_pn
contact_time_refinement_count
maximum_contact_time_refinement_depth
contact_residual_projection_count
maximum_contact_residual_projection_um
maximum_contact_complementarity_residual_pn_um
minimum_accepted_internal_wall_gap_um
accepted_negative_gap_count
contact_nonzero_velocity_zero_progress_count
contact_kinematic_interval_evaluations
contact_cumulative_position_path_um
contact_cumulative_velocity_path_um
contact_position_to_velocity_path_ratio
minimum_contact_interval_position_to_velocity_path_ratio
maximum_contact_interval_position_to_velocity_path_ratio
maximum_free_gap_kinematic_residual_um
directed_outlet_event_count
active_outside_lumen_violations
active_outside_accessible_domain_violations
opposite_wall_hydrodynamic_validity_warning_observations
```

For a valid new run, `accepted_negative_gap_count` and
`contact_nonzero_velocity_zero_progress_count` must both be zero, as must both
active-outside invariant counters. The complementarity residual must be small
at the configured residual tolerance, and the contact path ratio must approach one
under substep refinement. Time-refinement fields are administrative retry
audits; rejected positions, molecular state, and accepted-work extrema are not
merged into accepted results.

The metadata records `near_wall_xi_min`, `near_wall_xi_near`,
`near_wall_xi_far`, the radius-dependent minimum/maximum hydrodynamic
regularization gaps, and separate text definitions for the true and
hydrodynamic gaps. This makes it possible to verify that the mobility
regularizer did not become a geometric wall. The archived `run_config.yaml`
remains the authoritative record of all configured inputs; no separate
Numba-cache timing is exposed. Cardiac runs additionally use perfusion identifier
`deterministic_pulsatile_flux_halton_empty_start_v7` and store waveform, phase,
propagation, modulation-range, injection-rate-range, and cycle-count metadata.

Molecular runs additionally record target mode and density, state semantics,
maximum expected bond count, tension, force and torque, mean-field warnings,
and reaction-geometry diagnostics. Three limiter counters have different
meanings and must not be combined:
`molecular_bell_saturation_rhs_observations` and
`molecular_formation_rate_float_overflow_rhs_observations` count RHS
evaluations, whereas `molecular_capacity_limited_accepted_step_observations`
is incremented only when an accepted molecular state step is committed. The
console prints warnings for Bell saturation and accepted-step capacity
limiting; the formation-rate overflow count is retained in metadata without its
own warning.

The visualization/result loader has an explicitly read-only compatibility path
for historical Revised-v13 trajectories and still older files. It recognizes
`wall_contact_integrator=revised_v13_stateless_single_wall`,
`record_wall_slide_failure_duration_s`, and old `penetration_tolerance_um`
metadata so archived results can be inspected faithfully. The molecular pilot
may also label the v13 saved-RHS/realized-motion lock predicate when an archived
file is read. None of these branches is called by the current writer or v14
transport path, and their presence must not be interpreted as permission to
use v13 velocity projection, position repair, or negative accepted gaps in a
new run.

### 11.3 Revised-v8/v10 Candidate and Selected-Target Files

`molecular_target_candidates.npz` is written only by candidate-preparation mode.
It stores one `unit_id_grid` plus ragged offset/value arrays for unit segments,
unit children, and candidate memberships. This compact representation avoids
duplicating an `(nx, nz)` Boolean array for every candidate. It also stores the
physical `x_um`/`z_um` axes; lumen, solid-wall, open-boundary and candidate-
support masks; unit parent/root IDs, flow, length, volume, cylindrical
endothelial wall area, wall-area centroid/second moment, topology depth and
perfusion state; candidate IDs/kinds/parents/depths; per-root and network flow
fractions, volume, residence time, wall-area-weighted mean WSS, wall-area
fraction, centroid, radius of gyration, and expected bubble visits; the
finite-size injection rate and formal observation time; and the unresolved
junction-trace fallback count.

Schema v3 additionally stores `wall_area_weight_um2`,
`wall_segment_id_grid`, `accessible_wall_mask`,
`expected_bubble_visits_by_unit`, and the mapped/unmapped theoretical wall-area
totals. Schema v3 is required for Revised-v10 automatic generation. The loader
still accepts schema v1 and v2 for the manual candidate-tree workflow, but both
are deliberately rejected by the automatic workflow because they lack the
physical wall quadrature needed for area-weighted thresholding. The JSON
companion is an inspection index, while the NPZ remains authoritative.

The interactive selector writes `selected_molecular_target_mask.npz` by
default. Both workflows share the formal fields:

| Field                        |                Shape | Meaning                                                                                                  |
| ---------------------------- | -------------------: | -------------------------------------------------------------------------------------------------------- |
| `x_um`, `z_um`           | `(nx,)`, `(nz,)` | Strictly increasing physical cell-centre axes in micrometers.                                            |
| `target_mask`              |         `(nx, nz)` | Final Boolean target ROI; open caps are false. In automatic mode it contains wall-only positive patches. |
| `candidate_schema_version` |               scalar | Candidate-catalog schema that produced the selection.                                                    |
| `selection_mode`           |               scalar | `manual` or `automatic_spatial_heterogeneity`.                                                       |
| `target_mask_semantics`    |               scalar | Reminder that only eligible solid-wall intersections become molecular targets.                           |

Manual output additionally stores `selected_candidate_ids`, requested/achieved
wall-area values when applicable, and the historical manual-modification flag.
Its candidate union includes owned lumen and eligible wall support so the
selected vascular bed remains visually intelligible.

Automatic output instead stores `influence_region_mask`,
`influence_wall_mask`, the anchor candidate ID and physical centre/radius,
requested and achieved influence/positive wall-area fractions, positive
network wall-area fraction, correlation length, seed, mode count, field
algorithm, exact wavevectors/phases, threshold, sparse influence-wall field
values and flat indices, and the connected patch count. Its `target_mask` is
already wall-only. During a formal run, both output forms are intersected with
the current eligible solid wall, so lumen cells never become force-producing
molecular sites.

### 11.4 `molecular_target_field.npz`

This file exists only when `molecular_target.enabled=true` and
`output.save_npz=true`.

| Field                                                                     |                              Shape | Meaning                                                                                                                                     |
| ------------------------------------------------------------------------- | ---------------------------------: | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `enabled`                                                               |                             scalar | Saved target-field enable flag.                                                                                                             |
| `region_mode`                                                           |                             scalar | `mask_npz`.                                                                                                                               |
| `target_density_molecules_per_m2`                                       |                             scalar | Fixed target density assigned to selected wall sites.                                                                                       |
| `spacing_um`                                                            |                             scalar | CFD/particle grid spacing in micrometers.                                                                                                   |
| `x_coordinates_um`, `z_coordinates_um`                                |               `(nx,)`, `(nz,)` | World coordinates of the target-field grid cell centres.                                                                                    |
| `solid_wall_mask`                                                       |                       `(nx, nz)` | Grid overlay of eligible impermeable solid-wall sites; authoritative open faces are excluded.                                               |
| `target_wall_mask`                                                      |                       `(nx, nz)` | Target-positive subset of`solid_wall_mask`.                                                                                               |
| `target_density_field_molecules_per_m2`                                 |                       `(nx, nz)` | Target density on selected wall samples and zero elsewhere.                                                                                 |
| `open_boundary_mask`                                                    |                       `(nx, nz)` | Original inlet/outlet open-boundary mask retained for registration and invariant checks.                                                    |
| `wall_coordinates_xz_um`                                                |                         `(W, 2)` | Physical centres of the`W` authoritative face segments in \(\Gamma_w\); no artificial inward shift is applied.                            |
| `wall_normal_xz`                                                        |                         `(W, 2)` | Inward face normal for each row of`wall_coordinates_xz_um`.                                                                               |
| `wall_axis`                                                             |                           `(W,)` | Axis code of each authoritative solid face segment.                                                                                         |
| `wall_length_um`                                                        |                           `(W,)` | Physical segment length used by the analytic target-reaction interval integration.                                                          |
| `wall_target_positive`                                                  |                           `(W,)` | Boolean target flag mapped once from the wall-cell overlay to each authoritative solid face; it avoids runtime half-grid nearest-cell ties. |
| `source_mask_npz_path`                                                  |                             scalar | Present only for`mask_npz`; absolute resolved path of the coordinate-aware source ROI used by this run.                                   |
| `region_x_coordinates_um`, `region_z_coordinates_um`, `region_mask` | `(mx,)`, `(mz,)`, `(mx, mz)` | Present only for`mask_npz`; normalized copy of the coordinate-aware source ROI.                                                           |

Here `W` is the number of authoritative solid face segments, which need not
equal `count_nonzero(solid_wall_mask)` near corners. The face centres, axes,
lengths, normals, and per-face target flags are the molecular wall geometry; the
masks remain the grid-registration overlay used by selectors and visualizers.
The reaction-area routine projects the particle centre directly onto those finite segments
and integrates only target-positive face intervals, ensuring
\(\Gamma_T\subseteq\Gamma_w\) without erasing legitimate side-wall target
support near an opening.

The loader also accepts the source aliases `x_um`, `z_um`, and `target_mask`
when configured or found as fallbacks, but the saved output always uses the
normalized key names above.

### 11.5 Molecular Contact-Pilot YAML

Each `molecular_contact_pilot_*_capture_ratio_*.yaml` report records the report
kind `molecular_contact_exposure_and_predeclared_da_on_sweep`, interpretation,
unit conversion, capture/rest-length ratio, aggregate
pilot summary, per-bubble cumulative contact, contiguous contact events, and
the Cartesian product of the configured `Da_on`, target-density, and
ligand-density sensitivity levels. Its study context also records the target
source, target/wall sample counts, trajectory schema and record counts,
analytic area method, and invalid-normal observation count. A zero or invalid
wall normal cannot create contact area. The report also records realized-centre
velocity semantics, right-censored events, and any target-contact time that
overlaps the applicable numerical-lock predicate. For current v15 files this is
the last accepted internal constrained velocity versus same-ID realized
saved-frame motion; the v13 predicate is retained only when an archived v13
file is read. Every scenario row records the predeclared
`da_on_reference_time_s` and association-parameter source
`predeclared_dimensionless_scenario_independent_of_observed_exposure`. A
missing report can mean the sweep was disabled. Zero observed target contact is
a valid exposure result and does not suppress or alter the independently
predeclared scenario table.

## 12. Flow Diagnostics

`flow_diagnostics.yaml` contains:

```text
stage_divergence:
  div_initial
  div_after_boundary_apply
  div_after_projection
  div_final_saved

divergence_stats:
  core
  boundary
  lumen

mask_quality:
  n_components
  hole_cells
  narrow_cells
  narrow_fraction
  narrow_in_inlet_or_outlet
  narrow_in_junction_core
  min/p1/p5/median resolved diameter

outlet_flux:
  n_outlets
  max_relative_error

narrow_segments:
  narrow_segments_count
  narrow_segments_fraction
  segments_count
```

The divergence heatmap overlays final divergence with lumen contour, inlet,
outlet, wall, and segment ID. The outlet CSV is the authoritative per-outlet
flux check because it is computed from the final saved face fluxes.

## 13. Visualization

The implementation has three related visualization paths. They do not expose
the same interaction or data semantics and should not be treated as
interchangeable.

### 13.1 Automatically Generated PyVista/VTK CFD Artifacts

`generate_microbubble_trajectories.py` automatically renders an initial and a
final CFD scene after the accepted field and particle data have been saved.
Each scene is a linked `1 x 2` PyVista canvas:

| Panel | Background                                | Overlay                                                                                                        |
| ----- | ----------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| Left  | Velocity-magnitude field with LIC texture | Validated continuous root-inlet-to-terminal-outlet stream tubes, colored by speed                              |
| Right | Projection-pressure field                 | Spatially thinned velocity glyphs plus a translucent white copy of the representative validated streamline set |

The initial scene uses the initialized velocity and a zero pressure reference.
The final scene uses the accepted velocity and final projection potential. Axis
geometry is physical X-Z geometry and is not stretched independently. The
offline HTML files allow camera zoom, pan, and rotation, but they are exported
scenes rather than live numerical probes.

The formal `*_streamlines.vtp` files contain only paths that pass the endpoint,
direction, topology, and root-to-outlet continuity checks. Rejected VTK trace
candidates are kept in `*_streamlines_diagnostic.vtp`, and the corresponding
CSV is the human-readable continuity audit. A missing required outlet path is a
rendering failure, not permission to silently display an incomplete formal
streamline set.

The formal VTP retains **all** validated paths and their full integration-point
data. The automatically exported HTML/PNG scene uses a separate display-level
subset: by default it selects two representative paths per outlet and then adds
the shortest validated paths needed to retain vessel-topology coverage. Only
this display copy may be polyline-decimated (target reduction `0.60`, accepted
only when every simplified segment remains in fluid); the formal VTP is never
decimated. Therefore the number of tubes visible in an HTML scene is not the
number of validated paths saved in the formal streamline dataset.

The native `*_flow_field.vti` files use VTK `ImageData` with dimensions
`(nx + 1, nz + 1, 1)`. VTK X is physical X, VTK Y is physical Z, and VTK Z is
zero. Cell arrays are flattened in Fortran order and include the masks, vessel
ID, wall distance, speed, projection pressure, X-Z velocity, inlet/outlet
labels, wall-shear display mask, and stage-appropriate WSS. The point-data
`lic_intensity` array supplies the LIC texture. Initial WSS values are NaN;
initial pressure is the zero reference.

Formal and diagnostic VTP paths carry sampled velocity, speed, pressure,
integration time, wall distance, sampled vessel IDs, path IDs, path progress,
termination information, intended and reached outlets, and completeness flags.
The continuity CSV contains one row per traced path with these identifiers plus
point count, arc length, endpoint distances, direction-alignment fraction, and
sampled/covered vessel IDs.

The wall-shear CFD artifact is generated only for the accepted final flow. It
shows the closed solid-wall band and the existing nonnegative 2D in-plane WSS
proxy in Pa; it is not an exact three-dimensional CFD traction solution.

### 13.2 Live Numerical Probe with trame

The trame viewer reads the saved native VTI/VTP files and adds live cell-value
probing. Its defaults are `stage=final`, `view=flow`, host `127.0.0.1`, and port
`8080`. It does not open a browser unless `--open-browser` is supplied. Both
initial and final flow views are valid, but `wall-shear` is final-stage only.
For responsiveness, the live flow view rereads the full formal VTP and selects
up to two representative paths per outlet; unlike the automatic renderer's
display subset, this live-view selection does not add topology-coverage paths.

```powershell
D:\anaconda3\envs\pmp\python.exe -m ulm_microbubble_traj_gen.vis_utils.trame_flow_viewer --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --stage final --view flow --open-browser
```

```powershell
D:\anaconda3\envs\pmp\python.exe -m ulm_microbubble_traj_gen.vis_utils.trame_flow_viewer --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --stage final --view wall-shear --open-browser
```

### 13.3 Microbubble Trajectory Visualization

The trajectory entry point is:

```text
ulm_microbubble_traj_gen/visualize_microbubble_flow.py
```

It always regenerates `lumen_holes.png`, `narrow_lumen_cells.png`,
`initial_flow_field.png`, and `final_flow_field.png`. A normal invocation also
creates a GIF. `--no-animation` creates the default final-frame snapshot when
no explicit `--snapshot` path is provided. `--snapshot` always requires a path.
An MP4 is selected by giving an `.mp4` output path and requires an available
FFmpeg writer. If `--result-dir` is omitted, the loader selects the newest
complete result folder by directory modification time, not by parsing its name;
`--list` shows all folders containing both required NPZ files.

The four Matplotlib quality maps do not reuse the validated VTK streamlines.
Their initial/final flow panels call Matplotlib `streamplot` locally on a
subsampled grid with density `1.35`, so those decorative lines do not carry the
formal root-to-outlet continuity guarantee. Animation defaults to `--stride 3`
and explicitly appends the final saved frame when the stride misses it. Use
`--stride 1` when every stored simulation frame must be rendered.

By default, `--max-bubbles 0` draws the complete active per-frame population.
Positive values deliberately limit the display count. Tails follow permanent
bubble IDs, so records from different IDs are never joined even as the ragged
active population changes between frames.

| Color mode     | Bubble color means                                                                                                                                                                             |
| -------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `speed`      | Realized centre speed derived from same-permanent-ID accepted-position differences divided by saved`dt_s`; it does not use a replacement-ID jump or the instantaneous constrained RHS speed. |
| `wall_shear` | Sampled nonnegative wall-shear proxy in Pa.                                                                                                                                                    |
| `vessel_id`  | Saved 1-based nearest-vessel ID.                                                                                                                                                               |
| `diameter`   | Immutable physical bubble diameter in`um`.                                                                                                                                                   |
| `clearance`  | Finite-radius bubble-surface wall gap in`um`.                                                                                                                                                |
| `active`     | Compatibility Boolean. Current ragged continuous-perfusion records are all active and padded inactive lanes are filtered, so visible bubbles normally all have value`1`.                     |

The opaque bubble core is a physical-size disk drawn at the saved diameter.
The default three-layer translucent halo is a visual aid only, with
`--glow-scale 3.0`. Its apparent extent may cross a wall without implying that
the physical core penetrated it. `--no-glow` disables only the halo and does
not change particle size or transport. When full mobility diagnostics are
available, the small internal marker follows the saved `+Y` rotation angle.
Wall contact is shown by a separate orange ring outside the visibility halo;
for a current Revised-v15 file any negative gap uses red, so neither state
overwrites the scalar color of the physical core. Historical v13 and older
files use their archived read-only warning convention. The frame title reports
the active contact count.
The current trajectory renderer still does not offer color modes for angular
velocity, collision force, gap ratio, near-wall weight, or two-wall warnings.
When `molecular_target_field.npz` is present, it overlays the exact
target-positive wall cells plus a subtle visibility halo and density label. It
does not yet offer bond-state color modes; molecular trajectory arrays remain
numerical NPZ outputs for separate analysis.

## 14. Implementation Assumptions

| Assumption                                                                                        | Why the code needs it                                                                                                                                                                                                                                                                                                                  | Consequence                                                                                                                                                                                                                                            |
| ------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Exported`flow_rate` is in `um^3/s`.                                                           | Boundary flux and initial speed use this unit.                                                                                                                                                                                                                                                                                         | `q2d` is `um^2/s`; velocity is saved in `um/s`.                                                                                                                                                                                                  |
| Exported`radius` is in `um`.                                                                  | Geometry, mask diagnostics, initial profile, and particle checks use micrometers.                                                                                                                                                                                                                                                      | Radius and grid spacing must share units.                                                                                                                                                                                                              |
| The field is X-Z planar.                                                                          | The solver grid has dimensions`x` and `z`.                                                                                                                                                                                                                                                                                         | Y is fixed; this is not a full 3D flow volume.                                                                                                                                                                                                         |
| Parent-to-child vessel direction approximates local flow direction.                               | Initial velocity and open-boundary orientation use vessel direction.                                                                                                                                                                                                                                                                   | Flow direction quality depends on the exported vessel tree orientation.                                                                                                                                                                                |
| The raster lumen is the 2D domain\(\Omega\), while authoritative CFD faces classify its boundary. | Flow and particle lifecycle must share open faces instead of reconstructing different boundaries from labels.                                                                                                                                                                                                                          | Holes or disconnected components stop the run; open faces are excluded from\(\Gamma_w\), and directed outlet sections govern termination.                                                                                                              |
| Solid walls are impermeable face boundaries.                                                      | Wall flux is fixed to zero in the face-flux projection.                                                                                                                                                                                                                                                                                | No post-projection wall cleanup is performed.                                                                                                                                                                                                          |
| Open boundaries are mask boundary faces.                                                          | Inlet/outlet fluxes are prescribed on faces, not thick cell bands.                                                                                                                                                                                                                                                                     | Final flux diagnostics integrate the same face fluxes used by the solver.                                                                                                                                                                              |
| A mobility-mode microbubble is a finite-size rigid sphere at low Reynolds number.                 | The reduced state retains X-Z center translation and rotation about`+Y`; inertia is omitted.                                                                                                                                                                                                                                         | The default result is not a passive streamline, deformable bubble, or full 3D trajectory.                                                                                                                                                              |
| Particle-to-flow coupling is one way.                                                             | The steady field is solved once before particle transport.                                                                                                                                                                                                                                                                             | Particles, pair forces, and concentration-driven perfusion do not modify the fluid field.                                                                                                                                                              |
| The accepted steady field is the cardiac-cycle mean reference.                                    | The runtime modulates the accepted field without rerunning the pressure projection at every particle stage.                                                                                                                                                                                                                            | Cardiac output is a retarded-phase quasi-steady kinematic approximation, not a compliant transient CFD solution.                                                                                                                                       |
| The synthetic ECG envelope is only a positive flow surrogate.                                     | No measured subject-specific inlet waveform is currently supplied.                                                                                                                                                                                                                                                                     | BPM, phase, propagation speed, and surrogate choice must be reported; electrical ECG amplitude is not interpreted as physical pressure.                                                                                                                |
| The nearest wall is locally planar for hydrodynamic coefficients.                                 | The cited mobility is a single-sphere, single-planar-wall model.                                                                                                                                                                                                                                                                       | Curvature and opposing-wall hydrodynamic corrections are not modelled. Revised v15 rejects a true simultaneous distinct-wall contact instead of choosing a normal or adding an unsupported second constraint.                                          |
| A rigid closed wall supplies a centred, frictionless, one-sided normal reaction.                  | Revised v15 must prevent finite-size overlap without inventing tangential friction or a primary post-step geometric repair.                                                                                                                                                                                                            | Reaction is predicted in closed form through the complete local generalized mobility; molecular bonds remain separate loads, and every accepted state satisfies\(g_R\ge0\).                                                                            |
| The true gap and hydrodynamic regularization are different quantities.                            | Near-wall coefficient formulas require a finite algebraic lower scale at contact.                                                                                                                                                                                                                                                      | Only\(h_{\mathrm{hyd}}\) is regularized; geometry, molecular capture, outlet lifecycle, and saved clearance always use the unregularized \(g_R=d_w-R\).                                                                                                |
| A finite-size centre can use only the inlet-connected clearance component.                        | Local positive gap alone cannot prove that a radius can reach a branch from an inlet.                                                                                                                                                                                                                                                  | Every admitted and continuing centre must lie in\(\Omega_R^{\mathrm{in}}\); an inaccessible endpoint is a hard error.                                                                                                                                  |
| Pair hydrodynamics are block diagonal except for the short-range engineering collision force.     | No long-range many-body mobility tensor is assembled.                                                                                                                                                                                                                                                                                  | Distant particles have no hydrodynamic interaction and soft overlap can occur during collision relaxation.                                                                                                                                             |
| Local wall shear for particle loading comes from the signed 2D velocity gradient.                 | The mobility model needs direction as well as magnitude.                                                                                                                                                                                                                                                                               | The saved nonnegative WSS proxy is a trajectory diagnostic, not the force input.                                                                                                                                                                       |
| Collision relaxation time is a model parameter.                                                   | Pair stiffness is derived from relative mobility and this time scale.                                                                                                                                                                                                                                                                  | Keep it fixed when comparing`dt`, `dt/2`, and `dt/4`.                                                                                                                                                                                            |
| Continuous perfusion converts an imposed concentration into a planar injection rate.              | The rate is the inlet number concentration multiplied by accepted finite-size section flux.                                                                                                                                                                                                                                            | It is concentration-calibrated within the 2D/2.5D model, not a full three-dimensional bolus-injection model.                                                                                                                                           |
| A molecular target is a reviewed physical-wall ROI with fixed density.                            | Revised v8 provides topology/flow-based beds for manual selection. Revised v10 uses one accessible bed only as a coarse anchor, then thresholds a seeded spatially correlated field over physical-area-weighted wall sites inside a compact tissue-space region. Only the coordinate-aware Boolean mask is consumed by the formal run. | Automatic placement controls influence size, positive fraction, correlation scale and realization seed; it does not infer biology. Selection provenance and world-coordinate registration are reproducible inputs; density is not depleted by binding. |
| Molecular binding is a deterministic local mean field.                                            | The state stores expected bond count and total extension rather than individual random bonds.                                                                                                                                                                                                                                          | Low expected molecular counts weaken this approximation; the configured warning threshold is diagnostic only.                                                                                                                                          |
| One ligand-target pair and one local planar reaction disk represent the molecular interface.      | The current kinetics and mechanics use one density pair, one Hookean bond law, and nearest-wall geometry.                                                                                                                                                                                                                              | Multiple receptor species, shared target competition, catch bonds, curved surfaces, and opposing-wall reaction geometry are not represented.                                                                                                           |
| Molecular force is one-way coupled through the same particle mobility operator.                   | Bond and collision loads alter particle translation/rotation after the CFD field is accepted.                                                                                                                                                                                                                                          | Adhesion does not feed back into velocity, pressure, wall shear, or target density.                                                                                                                                                                    |

## 15. Common Reproduction Mistakes

| Mistake                                                                                                                                          | Correct handling                                                                                                                                                                                                                                                                                                                                                            |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Running the microbubble generator before exporting vessels.                                                                                      | First run`ulm_vascular_model_generator/vessel_generation.py`.                                                                                                                                                                                                                                                                                                             |
| Using only SWC for physics.                                                                                                                      | Use`.vessels.npz`; SWC does not carry all required flow/viscosity fields.                                                                                                                                                                                                                                                                                                 |
| Expecting the old legacy diffusion-only flow.                                                                                                    | The active solver is finite-volume face-flux projection with PhiFlow viscous relaxation and SciPy CG pressure correction.                                                                                                                                                                                                                                                   |
| Calling the solver`phiflow_incompressible_2d`.                                                                                                 | Use`phiflow_viscous_fv_projection_2d`; the misleading name is rejected.                                                                                                                                                                                                                                                                                                   |
| Expecting thick inlet/outlet velocity bands.                                                                                                     | Open boundaries are selected mask boundary faces.                                                                                                                                                                                                                                                                                                                           |
| Reconstructing particle openings from`inlet_label`/`outlet_label`.                                                                           | Current production uses the authoritative saved`open_face_*` and `open_section_*` geometry. Labels are diagnostics, not an alternative \(\Gamma_w\) or outlet definition.                                                                                                                                                                                               |
| Ignoring nonconvergence.                                                                                                                         | Nonconverged flow raises`FlowConvergenceError`; trajectories are not generated.                                                                                                                                                                                                                                                                                           |
| Expecting wall cleanup after projection.                                                                                                         | Wall no-penetration is enforced by face-flux constraints; no final velocity overwrite is allowed.                                                                                                                                                                                                                                                                           |
| Treating`narrow_cells` as near-wall cells.                                                                                                     | Narrow means`D_eff_px < 8`; near-wall behavior is a different concept.                                                                                                                                                                                                                                                                                                    |
| Interpreting quick-test output as a production result.                                                                                           | Quick-test uses coarse`8.0 um` spacing and reduced runtime controls; it is a plumbing test, not a wall-shear or mobility-quality run. Preserve the `--quick-test` flag because copied `run_config.yaml` does not contain the in-memory overrides.                                                                                                                     |
| Assuming a passive, fixed-population, or Taichi particle mode remains selectable.                                                                | These implementations have been removed; public configuration and dispatch support only mobility-based continuous perfusion.                                                                                                                                                                                                                                                |
| Treating`bubble_diameter_um=2.0` as the current fixed diameter.                                                                                | The fixed-diameter option has been removed. Every permanent ID receives one immutable diameter from the required configured uniform`1.5..2.5 um` **prior** through the deterministic low-discrepancy sequence. Because finite-size-accessible inlet flux depends on radius, the marginal diameter distribution of scheduled events need not remain exactly uniform. |
| Treating`n_steps` as the stored frame count.                                                                                                   | The output has`n_steps + 1` frames because frame 0 is included.                                                                                                                                                                                                                                                                                                           |
| Expecting a configured fixed number of bubbles in every frame.                                                                                   | The former`n_bubbles` option has been removed. Population emerges approximately as injection rate times residence time and varies with exits and inlet waiting.                                                                                                                                                                                                           |
| Treating a flat-record position as a reusable particle lane.                                                                                     | Use`record_bubble_id`. Permanent IDs and state slots are monotone and never reused; frame slices are ragged active-observation lists.                                                                                                                                                                                                                                     |
| Reading the trajectory NPZ as dense`(frames, bubbles, ...)` arrays.                                                                            | Use`frame_offsets` to slice flat `record_*` arrays and use the registry arrays for lifecycles.                                                                                                                                                                                                                                                                          |
| Interpreting continuous-perfusion`record_velocities_um_s` as either fluid velocity, free pre-contact velocity, or realized frame displacement. | It is the last accepted internal v15 constrained velocity when available, with a saved-state RHS fallback only for a newly admitted ID. Use same-ID`record_realized_velocities_um_s` for saved-frame motion and optional `record_fluid_velocities_um_s` for the carrier sample.                                                                                         |
| Treating the accepted field NPZ as an instantaneous cardiac phase.                                                                               | It is the cycle-mean reference field. Use the saved waveform, delay field, and per-record multiplier to interpret the particle-time carrier flow.                                                                                                                                                                                                                           |
| Converting pulse delay with samples per cardiac period.                                                                                          | Delay is`root_path_distance / pulse_propagation_velocity` in seconds. Waveform evaluation is continuous and periodic at internal-stage time.                                                                                                                                                                                                                              |
| Keeping injection events uniformly spaced while enabling pulsatility.                                                                            | The runtime inverts the integrated time-dependent number flux, preserving inlet concentration under a pulsatile volume flux.                                                                                                                                                                                                                                                |
| Treating`wall_contact_threshold_um` as a clearance barrier or constraint trigger.                                                              | It labels saved near-wall observations only. The physical constraint is\(g_R\ge0\); exact contact at \(g_R=0\) is intentionally permitted.                                                                                                                                                                                                                                  |
| Treating`contact_geometry_tolerance_um` as an allowed negative gap or a positive stand-off.                                                    | It only bounds the small exact-gradient endpoint residual projection. Every accepted endpoint and complete chord still require\(g_R\ge0\).                                                                                                                                                                                                                                  |
| Using the hydrodynamic regularization gap for wall feasibility or molecular capture.                                                             | Only mobility coefficients use\(h_{\mathrm{hyd}}\). Admission, contact, saved clearance, and target reaction geometry use the unregularized \(g_R=d_w-R\).                                                                                                                                                                                                                  |
| Assuming a locally feasible centre is reachable from the inlet.                                                                                  | Require membership in\(\Omega_R^{\mathrm{in}}\); a positive local gap does not prove that the finite-size centre can pass every upstream bottleneck.                                                                                                                                                                                                                        |
| Treating`contact_max_time_refinements` as a displacement-halving count.                                                                        | It bisects true physical time. Both half intervals are recomputed and together cover the original duration; exhaustion is an error.                                                                                                                                                                                                                                         |
| Assuming exact swept-wall checks prove time-step convergence.                                                                                    | Piecewise-bilinear path analysis prevents geometric tunnelling, but does not resolve mobility, collision, or bond transients. Perform fixed-parameter time-step convergence studies.                                                                                                                                                                                        |
| Using projection as the primary rigid-wall solver.                                                                                               | Revised v15 predicts the closed-form mobility reaction before movement. Only a small endpoint roundoff residual may be projected along the exact gradient; a materially infeasible chord is recomputed as two real half-time intervals or fails explicitly.                                                                                                                 |
| Advancing beyond an outlet and checking lifecycle on the next step.                                                                              | Directed\(\psi_k\) crossing has first priority, the centre is committed exactly on the outlet plane, and the permanent ID terminates immediately.                                                                                                                                                                                                                           |
| Interpreting the unilateral reaction as molecular adhesion or wall friction.                                                                     | It is the minimum nonnegative normal force needed for rigid nonpenetration. Bond force/torque remains separate, and tangential motion is retained through the full mobility matrix.                                                                                                                                                                                         |
| Treating`xi_*`, collision layer, or collision relaxation defaults as universal literature constants.                                           | They are explicit engineering parameters and must be reported and sensitivity-tested.                                                                                                                                                                                                                                                                                       |
| Assuming`neighbor_search=auto` always means a cell list.                                                                                       | `auto` deliberately keeps the faster deterministic all-pairs kernel below 1024 active particles and switches to the sparse cell list only at or above that implementation threshold.                                                                                                                                                                                      |
| Interpreting the opposite-wall warning as a corrected double-wall solution.                                                                      | It is only a hydrodynamic-validity heuristic. A true simultaneous contact with distinct solid walls is a separate hard error outside the current single-wall model.                                                                                                                                                                                                         |
| Reusing the saved WSS scalar as the mobility shear load.                                                                                         | The particle load uses the signed local velocity gradient; WSS is separately sampled for output.                                                                                                                                                                                                                                                                            |
| Assuming a blocked inlet event is force-placed or that waiting is strict FIFO.                                                                   | Production continuous perfusion never inserts an overlapping fallback. Blocked IDs remain upstream and are retried in deterministic ID order, but a later geometrically admissible event may enter while an earlier one remains blocked. Inspect waiting counts and admission times.                                                                                        |
| Assuming a quoted YAML value such as`"false"` is Boolean false.                                                                                | The current loader applies Python truth conversion to several flags; a nonempty quoted string is truthy. Use unquoted YAML`false`/`true`.                                                                                                                                                                                                                               |
| Assuming a timestamp directory name is globally unique.                                                                                          | The timestamp has one-second resolution and the directory is created with`exist_ok=True`; two launches resolving to the same output root and second can share/overwrite files. Serialize launches or use distinct `output.results_dir` values.                                                                                                                          |
| Defining a molecular target with matrix indices or an image lacking physical axes.                                                               | Use a Boolean NPZ mask with strictly increasing physical X/Z axes in micrometers. Prefer the selector's canonical`x_um`, `z_um`, and `target_mask` output.                                                                                                                                                                                                            |
| Treating an automatically generated synthetic target as a diagnosed tumour region.                                                               | Revised v10 uses topology/flow only for an accessible coarse anchor, then creates physical-space correlated wall patches from configured influence size, positive fraction, correlation length and seed. These controls do not establish patient-specific molecular expression. Review the mask and preserve its provenance before assigning biological meaning.            |
| Assuming every ROI pixel becomes adhesive.                                                                                                       | The ROI is intersected with canonical solid face segments and authoritative open faces are excluded. Verify both`target_wall_mask` and the saved wall face geometry in `molecular_target_field.npz`.                                                                                                                                                                    |
| Entering densities in`molecule/um^2` into fields named `*_molecules_per_m2`.                                                                 | Convert explicitly:`1 molecule/um^2 = 1e12 molecule/m^2`. Only the scenario-sweep level fields are intentionally expressed per square micrometer.                                                                                                                                                                                                                         |
| Copying a literature forward rate in`s^-1` into `association_rate_m2_per_molecule_s`.                                                        | The implemented association coefficient is an effective two-dimensional rate in`m^2/(molecule s)` and needs compatible calibration or a clearly labelled `Da_on` sensitivity derivation.                                                                                                                                                                                |
| Interpreting`record_bond_count_expected` as an integer/binary bound state.                                                                     | It is a continuous deterministic expectation. The runtime has no discrete single-bond identity or categorical trapped state.                                                                                                                                                                                                                                                |
| Interpreting`record_bond_formation_rate_bonds_s` as the configured `k_on`.                                                                   | This output stores the realized formation flux in expected bonds per second; the configured coefficient has units`m^2/(molecule s)`.                                                                                                                                                                                                                                      |
| Treating positive target overlap as proof of adhesion.                                                                                           | `record_target_reaction_area_um2` and overlap fraction describe geometric molecular contact; inspect expected bond count and molecular load to assess modeled binding.                                                                                                                                                                                                    |
| Expecting target geometry alone to change trajectories.                                                                                          | With`molecular_binding.enabled=false`, the target is passive and is used only for saved geometry or the no-bond contact pilot.                                                                                                                                                                                                                                            |
| Treating`mean_field_warning_count` as a physical formation threshold.                                                                          | It only warns that expected molecular populations may be too small for a deterministic mean field; it does not turn kinetics on or off.                                                                                                                                                                                                                                     |
| Replacing`da_on_reference_time_s` with the no-bond pilot's observed contact time.                                                              | Fix the reference time before transport. The pilot reports exposure only;`k_on = Da_on/(rho_T*t_ref)` is constructed from predeclared scenario inputs even when observed contact is zero.                                                                                                                                                                                 |
| Treating a contact-pilot YAML as a calibrated or already executed adhesion result.                                                               | It is a hypothesis-generation sensitivity report. Merge a selected, justified scenario into a full config and run a separate bound-particle simulation.                                                                                                                                                                                                                     |
| Refining only internal substeps when validating contact duration.                                                                                | Internal substeps test dynamics, while the pilot's rectangle-rule contact duration is limited by saved-frame`particles.dt_s`; refine both independently.                                                                                                                                                                                                                  |
| Expecting 3D particle motion.                                                                                                                    | The current output has fixed Y and X-Z motion.                                                                                                                                                                                                                                                                                                                              |
| Expecting MAT outputs.                                                                                                                           | The active field-based pipeline writes NPZ outputs.                                                                                                                                                                                                                                                                                                                         |

## 16. Not Implemented in the Current Runtime

| Feature                                                                                                              | Status                                                                                                                                                                                                                                                                            |
| -------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Automatic inference, segmentation, or biological calibration of a disease target region                              | Not implemented. Revised v10 can deterministically generate a spatially heterogeneous synthetic target with configured wall-area fractions, physical correlation length and seed, but this is not tumour detection and does not supply biological density or kinetic calibration. |
| Discrete stochastic bond identities or categorical bound/dwell/trapped particle states                               | Not implemented. Revised v7 stores continuous expected bond count and total tangential extension; dynamic slowing or arrest is not a separate state flag.                                                                                                                         |
| Persistent per-bond wall-anchor coordinates or receptor identities                                                   | Not implemented. Existing mean bonds are reinterpreted in the current nearest-wall tangent/normal basis after motion.                                                                                                                                                             |
| Target depletion, internalization, diffusion, shared occupancy, multiple ligand-target pairs, or catch-bond kinetics | Not implemented. The current target density is fixed and each bubble uses the same single-pair deterministic mean-field law.                                                                                                                                                      |
| Bubble deformation, shell mechanics, acoustic radius response, rupture, coalescence, or fragmentation                | Not implemented; mobility particles are rigid spheres.                                                                                                                                                                                                                            |
| Acoustic radiation, buoyancy, Brownian diffusion, lift, added mass, history force, and particle inertia              | Not implemented.                                                                                                                                                                                                                                                                  |
| Long-range many-body hydrodynamic coupling                                                                           | Not implemented; only the short-range mobility-scaled soft collision force couples particles.                                                                                                                                                                                     |
| Exact curved-wall, bifurcation, or double-wall hydrodynamic mobility                                                 | Not implemented; the nearest local planar wall is used for coefficients. A proximity heuristic warns about possible opposing-wall influence, while a true simultaneous distinct-wall contact fails explicitly.                                                                    |
| Exact hard-sphere no-overlap constraint                                                                              | Not implemented; the relaxation model permits transient physical or collision-layer overlap.                                                                                                                                                                                      |
| Tree neighbour search                                                                                                | Not implemented. Deterministic all-pairs and sparse cell-list searches are implemented; no KD-tree is used in the hot collision loop.                                                                                                                                             |
| General adaptive dynamics subcycling or local-error-controlled`dt`                                                 | Not implemented. Revised-v15 complete-chord feasibility can bisect a failed physical-time interval, but this is not a general mobility/collision/bond accuracy controller.                                                                                                        |
| GPU/Taichi implementation of the mobility model                                                                      | Not implemented; the obsolete Taichi particle path has been removed.                                                                                                                                                                                                              |
| Full 3D volumetric flow                                                                                              | Not implemented; current grid is X-Z planar.                                                                                                                                                                                                                                      |
| Pure PhiFlow incompressible solver                                                                                   | Not implemented; pressure projection is SciPy finite-volume CG.                                                                                                                                                                                                                   |
| Legacy diffusion-only pseudo-flow                                                                                    | Removed from the active runtime.                                                                                                                                                                                                                                                  |
| Thick inlet/outlet velocity override bands                                                                           | Removed from the active runtime.                                                                                                                                                                                                                                                  |
| Post-projection wall velocity cleanup                                                                                | Removed from the active runtime.                                                                                                                                                                                                                                                  |
| 1D-3D strong coupling                                                                                                | Not implemented.                                                                                                                                                                                                                                                                  |
| Exact 3D no-slip wall reconstruction                                                                                 | Not implemented.                                                                                                                                                                                                                                                                  |
| Exact CFD wall shear stress                                                                                          | Not implemented; WSS is an in-plane proxy.                                                                                                                                                                                                                                        |
| Compliant-wall transient pulsatile CFD, pressure-wave reflection, and Womersley reconstruction                       | Not implemented. Cardiac transport uses a mean-field, retarded-phase kinematic modulation.                                                                                                                                                                                        |
| Subject-specific measured flow or pressure waveform                                                                  | Not implemented; the enabled default is a clearly labelled positive synthetic ECG-shaped surrogate.                                                                                                                                                                               |
| An inlet-concentration/number-flux pulsation control independent of carrier-flow pulsatility                         | Not implemented. With cardiac pulsatility enabled, event timing follows the same attenuated carrier-flow multiplier so the configured inlet concentration is preserved.                                                                                                           |
| Built-in trajectory overlays or color modes for molecular target and bond state                                      | Not implemented. Molecular target and bond arrays are saved numerically but are not consumed by the current Matplotlib trajectory renderer.                                                                                                                                       |
| MAT output for field-based trajectories                                                                              | Not implemented.                                                                                                                                                                                                                                                                  |

## 17. Reproduction Commands

Generate or refresh the vascular model:

```powershell
d:\anaconda3\envs\pmp\python.exe ulm_vascular_model_generator\vessel_generation.py
```

Run the full field-based microbubble generator:

```powershell
d:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\generate_microbubble_trajectories.py
```

Run a quick test:

```powershell
d:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\generate_microbubble_trajectories.py --quick-test
```

Run with a non-default config:

```powershell
d:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\generate_microbubble_trajectories.py --config ulm_microbubble_traj_gen\configs\physics_flow_config.yaml
```

Reuse a previously accepted field for repeated particle/parameter runs. The
source can be the result directory or the field NPZ itself:

```powershell
D:\anaconda3\envs\pmp\python.exe `
  ulm_microbubble_traj_gen\generate_microbubble_trajectories.py `
  --config ulm_microbubble_traj_gen\configs\molecular_contact_fixed_target_20s.yaml `
  --reuse-field-from ulm_microbubble_traj_gen\results\20260721_103516 `
  --skip-render
```

Prepare Revised-v8/v10 candidate vessel beds after solving and accepting the CFD
field. This mode does not generate particle trajectories:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\generate_microbubble_trajectories.py --config ulm_microbubble_traj_gen\configs\physics_flow_config.yaml --prepare-target-candidates
```

To make the same preparation run also write a deterministic, spatially
heterogeneous synthetic target, configure its influence size, positive-wall
fraction, physical correlation length, and reproducible realization explicitly:

```yaml
molecular_target_selection:
  default_mode: automatic
  influence_region_endothelial_wall_area_fraction: 0.10
  target_positive_wall_fraction_within_influence: 0.50
  target_correlation_length_um: 20.0
  random_seed: 42
  random_field_modes: 512
```

This writes `selected_molecular_target_mask.npz` plus
`automatic_molecular_target_selection.json`. The interactive selector can
rebuild and preview the same automatic realization. Candidate-tree checkbox
edits switch explicitly to the separate manual workflow rather than silently
modifying an automatic result.

Open the interactive hierarchy and export the selected Boolean mask:

```powershell
D:\anaconda3\envs\pmp\python.exe -m ulm_microbubble_traj_gen.vis_utils.trame_molecular_target_selector --result-dir "ulm_microbubble_traj_gen\results\<candidate-timestamp>" --open-browser
```

The default export is
`<candidate-timestamp>/selected_molecular_target_mask.npz`. Reference that file
from the complete formal-run YAML:

```yaml
molecular_target:
  enabled: true
  region_mode: mask_npz
  mask_npz_path: ../results/<candidate-timestamp>/selected_molecular_target_mask.npz
  mask_array_key: target_mask
  x_coordinates_key: x_um
  z_coordinates_key: z_um
```

The relative mask path is resolved against the directory containing the YAML.
Candidate preparation with `--quick-test` is useful only as a software smoke
test: its coarser/shorter CFD settings must not be used to report scientific
candidate metrics.

For a molecular run, use a **complete** YAML that contains the normal
input/domain/flow/particle sections, the selected `mask_npz` target, and the
justified molecular parameters. Files under
`configs/molecular_binding_scenarios/` are overlays and provenance records, not
stand-alone generator configs. When an effective two-dimensional association
coefficient has not been independently calibrated, declare the dimensionless
sweep axes and `da_on_reference_time_s` before transport. The empty-lumen
no-bond pilot then measures exposure; it does not calibrate or modify the
scenario table.

Run the focused hybrid geometry, contact, and lifecycle regressions:

```powershell
D:\anaconda3\envs\pmp\python.exe -m unittest ulm_microbubble_traj_gen.test_files.test_continuous_vessel_geometry_v16 ulm_microbubble_traj_gen.test_files.test_particle_predictive_contact_v15 ulm_microbubble_traj_gen.test_files.test_particle_hydrodynamic_fields ulm_microbubble_traj_gen.test_files.test_particle_inlet_flux_v14 ulm_microbubble_traj_gen.test_files.test_particle_constrained_step ulm_microbubble_traj_gen.test_files.test_particle_perfusion_time_transactions ulm_microbubble_traj_gen.test_files.test_contact_diagnostics_schema ulm_microbubble_traj_gen.test_files.test_particle_continuous_perfusion ulm_microbubble_traj_gen.test_files.test_molecular_transport_integration ulm_microbubble_traj_gen.test_files.test_molecular_pilot_runner -v
```

These tests cover authoritative open/solid face partitioning, directed outlet
planes, \(\Omega_R^{\mathrm{in}}\), exact finite-face true-gap consistency,
closed-form predictive complementarity, exact swept-disc forbidden pockets, outlet/wall tie
priority, exact outlet position/time, outside-domain failure, true physical-time
bisection, rejected-trial state isolation, contact and gap kinematics,
Python/Numba equivalence, current schema fields, historical read-only loading,
and v15 numerical-lock reporting.

Run the focused Revised-v7/v8/v10 regression tests:

```powershell
D:\anaconda3\envs\pmp\python.exe -m unittest ulm_microbubble_traj_gen.test_files.test_molecular_config ulm_microbubble_traj_gen.test_files.test_molecular_target_field ulm_microbubble_traj_gen.test_files.test_molecular_target_candidates ulm_microbubble_traj_gen.test_files.test_molecular_target_auto_selection ulm_microbubble_traj_gen.test_files.test_molecular_binding ulm_microbubble_traj_gen.test_files.test_molecular_contact_pilot ulm_microbubble_traj_gen.test_files.test_molecular_pilot_runner ulm_microbubble_traj_gen.test_files.test_molecular_scenario_configs ulm_microbubble_traj_gen.test_files.test_molecular_transport_integration -v
```

Record the exact pass count and package versions with each scientific run; the
test inventory can change as regression coverage is extended.

List available result folders:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\visualize_microbubble_flow.py --list
```

Create static quality maps and a final-frame snapshot for a specific result:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\visualize_microbubble_flow.py --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --no-animation --no-progress
```

Create a GIF with every active bubble visible in every stored frame:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\visualize_microbubble_flow.py --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --max-bubbles 0 --stride 1
```

Create an MP4 colored by physical diameter:

```powershell
D:\anaconda3\envs\pmp\python.exe ulm_microbubble_traj_gen\visualize_microbubble_flow.py --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --output "ulm_microbubble_traj_gen\results\<timestamp>\microbubble_flow.mp4" --color-mode diameter
```

Open the live final flow probe:

```powershell
D:\anaconda3\envs\pmp\python.exe -m ulm_microbubble_traj_gen.vis_utils.trame_flow_viewer --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --stage final --view flow --open-browser
```

Open the live initial flow probe:

```powershell
D:\anaconda3\envs\pmp\python.exe -m ulm_microbubble_traj_gen.vis_utils.trame_flow_viewer --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --stage initial --view flow --open-browser
```

Open the live final wall-shear probe (`initial wall-shear` is intentionally
unavailable):

```powershell
D:\anaconda3\envs\pmp\python.exe -m ulm_microbubble_traj_gen.vis_utils.trame_flow_viewer --result-dir "ulm_microbubble_traj_gen\results\<timestamp>" --stage final --view wall-shear --open-browser
```

Run validation tests:

```powershell
D:\anaconda3\envs\pmp\python.exe -m unittest discover -s ulm_microbubble_traj_gen\test_files -v
```

Record the exact discovered/pass count in the run log. Expected test output may
include deliberately triggered capacity, inlet-overlap, displacement-ratio,
contact-refinement, and VTK deprecation warnings.

Do not add `-p "test_*.py"` to that discovery command in the current environment:
trame may parse `-p` as its port argument when imported by visualization tests.
