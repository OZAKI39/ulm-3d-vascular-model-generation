# Reproduction boundaries

Root: /home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance
Python: /home/lzy/projects/computation_examples/_shared/python_environment/bin/python
The original external project is read-only. source_assets holds verified server snapshots.
Y_CLOT_REQUEST.json parameterizes position, footprint, height and physical constants. SOURCE_CONFIG.json gives external paths and identity transform. One request was executed.

Geometry entry: `python -m pd_clot.lammps_fragmentation.y_scene_task geometry --configuration verification/y_bifurcation_clot_placement_visual_integration_v1/Y_CLOT_REQUEST.json`.
Initialization entry: `python -m pd_clot.lammps_fragmentation.y_scene_task run0`; completed run0 is guarded against accidental repetition. Use a fresh output task directory for a separately authorized rerun; preserve this evidence.
Offline figure entry: `PYTHONPATH=verification/y_bifurcation_clot_placement_visual_integration_v1/visual_dependencies:. python -m pd_clot.lammps_fragmentation.y_scene_visualization`.
The figure entry reads saved files, traces the stored vector field only and never launches CFD or dynamics.
Installed export dependencies are isolated in visual_dependencies; shared / parallel environment unchanged.
The lightweight ZIP excludes native source VTU/VTP copies, full original background_flow.vtu and superseded draft renders; exact source paths/hashes remain in the register. Full local ParaView output is retained.
