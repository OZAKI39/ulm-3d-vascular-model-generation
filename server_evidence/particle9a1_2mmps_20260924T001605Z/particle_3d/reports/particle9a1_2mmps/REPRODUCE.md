# Reproduce this stage

Use the source identity recorded in `outputs/particle9a1_2mmps/provenance/*_config.json`.
All dynamics run on the server; the following post-processing is read-only and runs in WSL.

```bash
cd /home/lzy/projects/ulm_particle_3d_particle0
export PYTHONPATH=particle_3d/src PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1
.venv/bin/python particle_3d/scripts/analyze_particle9a1.py
.venv/bin/python particle_3d/scripts/audit_particle9a1_production.py
.venv/bin/python particle_3d/scripts/render_particle9a1.py
.venv/bin/python particle_3d/scripts/finalize_particle9a1.py
```

Server runner: `particle_3d/scripts/run_particle9a1.py`.
Stages must be `same12`, then gated `smoke30`, then gated `production`.
Inputs and source/event hashes must match before any stage cache is reused.
Original P65/P9-A comparisons are read from the protected previous diagnosis.
No new random admission events are sampled in this repair stage.
The new stage ledger retains all 500 current-flow admission events exactly.
The diagnostics observer defaults off. Same-12 and smoke use `--diagnostics`;
production retains velocity/omega/quaternion/gap/state for every accepted step,
without the additional full solver history overhead for each of 500 paths.
