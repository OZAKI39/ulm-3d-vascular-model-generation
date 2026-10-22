# Execution and reproduction

Local root: `/home/lzy/projects/computation_examples/_archives/mevo_observable_f10af_20261019/computation_examples/mevo_acoustic_microbubbles`. Branch: `dev/f10af-observable-multifidelity-20261019`. Server root: `/workspace/mevo_acoustic_microbubbles/f10af_observable_20261019`. Independent Python environments, immutable Git archives under execution/<commit>_<label>, unique Supervisor labels and bounded timeout. Actual execution commits are in `evidence/*_EXECUTION.json`, never replaced by the final documentation commit.

Activate only this worktree/snapshot with `source scripts/straight_tube_exp1/f10af/activate.sh`. This sets the separate AF environment and explicit current source import paths.

Read README before re-running: frozen plans and IDs must not be overwritten. Preserve original response and geometry files matching PREREGISTRATION protected hashes. Run `PYTHONPATH=src:. .venv/bin/python -m pytest tests/f10aa tests/f10ab tests/f10ae tests/f10af -q` on the configured server snapshot. Reproduce new experiments into a new output root, using protocol values and documented upstream files. The existing prepare/build functions intentionally reject already-created outputs. Full raw server files remain authoritative; exact accepted-frame datasets are provided locally for inspection.

## Evidence

- `PREREGISTRATION.json`
- `evidence/`
- `INPUT_OUTPUT_SHA256.json`
- `evidence/remote_environment.lock`

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false
