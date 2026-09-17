# 2D microbubble trajectory generation — code and reports

This snapshot contains the local microbubble generation and visualization code,
configuration, existing tests, research reports, and small static report figures.
It intentionally excludes animations, simulation datasets, cached flow fields,
exported vessel data, reference-paper PDFs, and local execution caches.

The 261 original files (12,780,546 bytes, approximately 12.2 MiB) are preserved
byte for byte. [SYNC_MANIFEST.csv](SYNC_MANIFEST.csv) records their paths, sizes,
and SHA-256 hashes. See [SYNC_REPORT.md](SYNC_REPORT.md) for validation details.

## Included files

| Path | Contents |
| --- | --- |
| `generate_microbubble_trajectories.py` | Trajectory generation entry point |
| `visualize_microbubble_flow.py` | Saved-result visualization entry point |
| `utils/` | Flow, geometry, particle transport, binding, and visualization code |
| `configs/` | Original simulation configurations |
| `test_files/` | Existing unit and regression tests |
| `check_files/` | Specifications, review summaries, revisions, and research reports |
| `figure_draw/` | Figure scripts and their documentation |
| `figure_draw/outputs/` | 10 PDF and 11 PNG static report figures; no animations or arrays |
| `../ulm_vascular_model_generator/` | Companion generator code and configuration required by imports |
| `../report_sp.md` | Collision, molecular adhesion, and RBC drift explanation |
| `../transport.md` | Microbubble mobility and transport explanation |

The companion generator includes its `utils/`, `configs/`, package initializer,
two entry points, and README. Its separate Git history, datasets, reference
papers, and unrelated development outputs are not imported.

## Environment and entry points

Run these commands from the repository root. The existing environment and
requirements files are copied as supplied, not newly validated lockfiles.

```bash
conda env create -f ulm_microbubble_traj_gen_2D/environment-dolfinx.yml
conda activate ulm-dolfinx
python -m pip install -r ulm_microbubble_traj_gen_2D/requirements-cfd-visualization.txt
python -m pip install -r ulm_microbubble_traj_gen_2D/requirements-acceleration.txt
python ulm_microbubble_traj_gen_2D/generate_microbubble_trajectories.py --help
python ulm_microbubble_traj_gen_2D/visualize_microbubble_flow.py --help
```

See [DOLFINX_BACKEND.md](DOLFINX_BACKEND.md) for the original backend notes.
No Git LFS download is needed for the newly synchronized files.

## Inputs required for simulation

Data is not part of this synchronization. Before running an experiment, generate
or separately supply the matching `.swc` and `.vessels.npz` vessel inputs and any
molecular target mask required by the chosen configuration. Visualization also
requires separately supplied simulation results. The placeholder
[`vessel_swc_models/README.md`](../ulm_vascular_model_generator/vessel_swc_models/README.md)
explains where vessel inputs are expected.

Original code, configurations, and reports retain historical paths and claims:

- Some documentation and figure scripts refer to the former package name
  `ulm_microbubble_traj_gen`. The synchronized package has the `_2D` suffix.
- `configs/physics_flow_config.yaml` contains an absolute molecular target path
  under the old package name. Copy the configuration and point
  `molecular_target.mask_npz_path` to an available input before a new run.
- `configs/typical_trajectory_reference_steady.yaml` names a historical vessel
  run. Supply that matching input to reproduce the corresponding experiment;
  substituting a different model is a different experiment.
- Earlier reports describe earlier numerical backends. Consult the current
  code and configuration for the behavior of a new simulation.

The synchronization does not rewrite numerical code, regenerate outputs, or
claim that full CFD experiments have been rerun.
