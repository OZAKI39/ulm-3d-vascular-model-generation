# Stage 5A-R2 evidence and reproduction

The reviewed starting commit is `7c39828d93509fff5551b7a16537ca3aa4ee9c23`.
The executed V6 numerical kernel is frozen at `07958c962`; V3 campaign source
is `5198d3c705c6a3735dcf566f53810274c40a6624`, and V4A entry/controller source
is `7f83f802c6cc435918ea25d69dcbb0f4f4dacce9`. A later archive commit adds
analysis and delivery only. Reused V3 runs keep their actual original identities.

Start with [the final report](../../results/stage5/motion_solver_r2/final/REPORT.md),
[separated gates](../../results/stage5/motion_solver_r2/final/MOTION_SOLVER_R2_GATE.json)
and [figure index](../../results/stage5/motion_solver_r2/figures/index.md).
The raw accepted/rejected states, immutable failure matrices, fresh JUnit and
sealed protocols are separate evidence layers. Figure smoothness is not a
numerical or physical validation criterion.

## Environment and inputs

The isolated executed project is `/workspace/mevo_acoustic_microbubbles_stage5a_r2`.
It uses Python 3.12.3 and `requirements-stage5a-r2.lock.txt`; install into a new
venv, not the old Stage4L/Stage5A environment. `python.sh` fixes BLAS/OpenMP
threads to one. Recorded CUDA wheels require a compatible host driver; do not
change the driver to match a lock without a separate infrastructure decision.
Package licenses, installation logs, CPU quota/affinity and device identity are
under `results/stage5/motion_solver_r2/environment/`.

Large input fields and some trajectory exports are external to Git. Retrieve
the exact files named in each run's `IDENTITY.json` and the final
`ARTIFACT_MANIFEST.json`; verify SHA256 before execution. The old Stage4L source
and original 3 s input snapshots are required at their recorded paths. The
standalone repository checkout is not a self-contained copy of those datasets.
The formal validator rejects changed source or mismatched protocol identities.

```bash
# Run inside a separate reproduction checkout with all recorded inputs present.
uv venv --python /usr/bin/python3 .venv
uv pip install --python .venv/bin/python -r requirements-stage5a-r2.lock.txt
# The project is imported by PYTHONPATH; the lock does not need an editable path.
```

Check the existing sealed protocol without rerunning the immutable freeze writer:

```bash
bash scripts/stage5a_r2/python.sh -c 'from scripts.stage5a_r2.validation_v4a import verify; verify(); print("Sealed source/test/environment identities match")'
```

## Experiments

Exact executed commands, arguments, measured wall times and per-run outcomes
are stored in `final/CAMPAIGN_V4A/RESULT.json`, with earlier V1/V3/V4 evidence
retained. Use those commands only in an isolated output copy. The original
campaign directory is immutable and the controller intentionally refuses to
overwrite it. Replaying a completed invocation creates a new timestamped run;
it does not extend or resume an earlier user-stopped campaign.

The selected V4A configuration uses maximum macro dt 5 us and time-error factor
0.25. Independent refinement uses 2.5 us and 0.0625. Both retain original KKT,
resistance and contact tolerances. `run_v4a.py` checks the whole formal contract
before entering the unchanged V6 loop. The 1 ms gate controls the 5 ms stage;
neither an installation nor a successful first step authorizes a long campaign.

Fresh test command/file hashes and JUnit paths are in `final/TEST_GATE_V4A.json`.
V4A reuses the twelve V4 holdout probes explicitly, because only entry routing
changed. Two formerly failing WALL states are separate development regressions,
not fresh holdouts. No additional tuning follows that freeze.

## Postprocessing

After all isolated timing processes finish, rebuild the original endpoint
equations and complete the audits:

```bash
bash scripts/stage5a_r2/python.sh -m scripts.stage5a_r2.endpoint_audit
bash scripts/stage5a_r2/python.sh -m scripts.stage5a_r2.manufactured_final
bash scripts/stage5a_r2/python.sh -m scripts.stage5a_r2.partial_comparisons
bash scripts/stage5a_r2/python.sh -m scripts.stage5a_r2.artifacts --server
```

The two independent scientific audits may run concurrently (three CPU workers
each, six in total, BLAS=1) after formal timings end. The endpoint audit includes
the last accepted state of incomplete/failed prefixes and labels that scope;
passing this local equation check never changes a failed full-duration status.

The figures/video were rendered in the recorded local environment; its separate
identity is `environment/RENDER_ENVIRONMENT.json`. With the exact HDF5 files
available locally, from the project root:

```bash
PYTHONPATH=src:. python -m scripts.stage5a_r2.diagnostic_figures
PYTHONPATH=src:. python -m scripts.stage5a_r2.validation_figures
PYTHONPATH=src:. python -m scripts.stage5a_r2.motion_visuals --video --duration 0.001
PYTHONPATH=src:. python -m scripts.stage5a_r2.video_audit
PYTHONPATH=src:. python -m scripts.stage5a_r2.delivery
PYTHONPATH=src:. python -m scripts.stage5a_r2.artifacts
```

The MP4 uses a uniform 60 fps presentation clock and last actual accepted states,
with sample age shown. It never interpolates coordinates. Its exact renderer
snapshot is retained separately when later static-label edits are made.
Raw ledgers are not altered by postprocessing; lossless gzip copies support
bounded Git delivery. The manifest records actual local/server paths and hashes
for excluded large assets, not fictitious download links.

The new backend is opt-in `constrained_lubrication_r2`. Continuous dense motion
uses the bounded S1/S0 V6 driver with endpoint implicit resistance and explicit
external forces; it is not a fully implicit acoustic-force integrator. Promotion
requires both 5 ms numerical/scientific gates and comparable complete-run 5x
speedup. Passing a numerical gate does not validate dense secondary-force physics.
