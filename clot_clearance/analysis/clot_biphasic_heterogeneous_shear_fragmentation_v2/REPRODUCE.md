# Reproduction

Project: /workspace/y_clot/project. Activate /workspace/y_clot/env.sh.

Sources: analysis/clot_biphasic_heterogeneous_shear_fragmentation_v2/. Runtime/source/input hashes: verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/SOURCE_FREEZE.json.
The final isolated library is under verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z/build. The capacity-only prior comparison library is under baseline_capacity. Build commands and upstream object hashes are retained.

Completed order: generate_biphasic_material.py → preflight_figures.py (P0 displayed) → build_capacity_baseline.py → build.py → uniform_probe.py old_uniform/new_uniform (separate libraries) → checks.py → path_audit.py → freeze_sources.py → run_cases.py --name uniform → run_cases.py --name biphasic → analyze.py → figures.py → render_videos.py → report.py → final_delivery.py.

Do not rerun formal cases in this output directory. The sole two permitted formal attempts are preserved. Reproduction that changes physical trajectories requires a new RUN and new authorization/scope; read-only plotting does not.

Read-only plotting on the existing server environment:
```bash
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /workspace/y_clot/project/analysis/clot_biphasic_heterogeneous_shear_fragmentation_v2
python analyze.py
python stage_comparison.py --key peak_load
python stage_comparison.py --time-ms 2.55 --nearest
python stage_comparison.py --key upper_drive_off
python preflight_figures.py
python figures.py
python render_videos.py --smoke
python render_videos.py
python report.py
```

Y_CLOT_PROJECT_ROOT overrides the root for replay in an environment containing numpy/scipy/pandas/matplotlib/PyVista/VTK/Pillow/FFmpeg. Native execution additionally needs mpi4py, the LAMMPS Python binding and hashed native build dependencies; LD_LIBRARY_PATH must select the isolated library. No whole-history dependency directory was duplicated.

Frame indices: videos/*_FRAME_TIMES.csv. Full FP64 snapshots: states/<case>/keys/*.npz. Every accepted scalar/graph and exact damage increments remain in each case. All geometry uses SI internally. Original IDs start at 1.
