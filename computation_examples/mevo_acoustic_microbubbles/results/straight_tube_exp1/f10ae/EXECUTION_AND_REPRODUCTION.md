# Execution and reproduction

Local branch: `dev/f10ae-targeted-shape-recovery-20261019`. Server root: `/workspace/mevo_acoustic_microbubbles/f10ae_target_20261019`. Isolated copied environment plus pinned optimizer installation; full freeze and logs in evidence/server_logs. Source is exported from clean Git commits into immutable `execution/<commit>_<label>` snapshots. Jobs use independent `mevo_f10ae_*` Supervisor labels and timeouts.

From a fresh isolated execution with read-only historical dependencies: `PYTHONPATH=src:. .venv/bin/python -m scripts.straight_tube_exp1.f10ae.prepare`, then `low_pressure`, `workflow`, and `post`. Existing IDs are immutable and hashed; create a fresh result root for a new campaign. `SOURCE_F10AE_COMMIT.txt` and per-run F10AE_INPUT.json provide the actual code/input identity.

Do not overwrite historical stage results or launch hardware. Portable datasets store exact saved accepted frames, with each file recording its sampling; original all-step HDF5 remains on the server. Videos map to actual raw indices. No fabricated frames or state interpolation.

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
