# F6 execution and analysis

Start from the project directory on branch
`dev/fast-pzt-f6-target-navigation-20261010`. The fixed F5 snapshot is
`c8ecf9deece1bb0e86b0d0f1c574b984a30ebbd0`. F6 adds modules; it does not edit F1–F5.
See `docs/FAST_PZT_F6_DEVELOPMENT.md` for the actual delivered results and limits.

## Actual compute layout

`remote.py` targets SSH alias `mevo-new`, the user's new host at port15283.
It commits no source automatically and refuses dirty executable F6 files. Each
job extracts `git archive <HEAD>` to the isolated server folder
`/workspace/mevo_acoustic_microbubbles/f6/.execution/<commit>` and uses the existing
`/workspace/mevo_acoustic_microbubbles/.venv` (Warp1.18.0+cu12). Results are shared
under `f6/results/fast_pzt/f6`; immutable source markers and SHA-256 are in IDENTITY.
The older F5 remote launcher still targets the retired host and is not used.

The main inputs are the original F5 `inputs/` and the exact archived F2 mature
and F5 switch trajectories. The sync baseline includes them. F1–F5 server result
links are read-only inputs in this workflow, not a reason to rerun old simulations.

## Reproduction commands

Choose **new** labels; existing attempts and supervisor configurations are retained.
These commands launch real jobs, so they should only be used when a new execution
is intended. GPU admission and science must be serialized on the single4090.

```bash
# Run from computation_examples/mevo_acoustic_microbubbles; source already committed.
python3 -m scripts.fast_pzt_f6.remote admission --label admission_repro --timeout 300
# Wait for supervisor mevo_fast_pzt_f6_admission_repro to finish and inspect its checks.
python3 -m scripts.fast_pzt_f6.remote science --label navigation_repro --timeout 2400 \
  --checks admission_repro_checks
```

The science command reads the tracked `configs/fast_pzt/f6/protocol.json`; it never
reselects targets. Its legacy `--protocol` label is not a configuration override.
`preflight.py` records static target/field probes if a new protocol is deliberately
developed; it is not required to reproduce the frozen campaign. The archived
preflight draft predates final feedback, timestep and leakage-safe group additions;
`navigation_v1/REGISTERED_PROTOCOL.json` is the actual formal protocol.

The two independent CPU radial references may run while the GPU is occupied:

```bash
python3 -m scripts.fast_pzt_f6.remote reference --label reference_repro --timeout 300 \
  --protocol preflight_v2
```

This pins official APECSSv1.8 and records its dependency/build/modified example
driver hashes. Only the excitation callbacks change; radial library code is intact.
No CFD, acoustic basis or production response table is regenerated.

After the science process has closed all HDF5 files, archive to WSL with rsync
excluding `dependencies/` (third-party build products), then verify the manifest.
Do not copy open HDF5 files as final evidence. Local analysis needs NumPy, SciPy,
h5py and Matplotlib; Warp is not needed:

```bash
PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1 python -m scripts.fast_pzt_f6.analyze --label navigation_v1
PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1 python -m scripts.fast_pzt_f6.reference_audit --label reference_v2
```

Both analysis commands refuse to overwrite their output directories. They do not
alter trajectories. Repeating them requires a separate copy/output folder, not
deleting the first report. Admission scripts are bounded checks; no full old
campaign, random scan, stress test or decorative animation is needed.

## Data access and interpretation

```python
from mevo_mb.fast_sim.ai_samples import read_group_sample, validate_partitions
sample = read_group_sample("results/fast_pzt/f6/navigation_v1/zero_flow_feedback/trajectory.h5",
                           frame=100, history_frames=3)
```

HDF5 retains the inherited `positions_m`, `time_s`, `metrics`, particle IDs and
configuration alongside `source`, `medium`, `flow`, `bubbles`, `control`,
`trajectory`, `morphology`, `quality`. Hard links avoid copying inherited datasets.
The reader returns a short history of positions, group tensors/Rg/labels, local
flow velocity/gradient, wall gap, activity and actual control events, plus fixed
radii/IDs and source/quality metadata. Off-wall NaN wall distance means no active
wall model, not a measured clearance. All quantities are SI unless named otherwise.
`environment` is the authoritative scenario configuration. The inherited summary
has generic F5 flow/wall capability labels even when those features are disabled;
use the actual flags/provider, not those generic labels. `applied_velocity_m_s`
is explicitly zero when flow is disabled, while raw diagnostic NaNs and validity
masks are preserved. Run-level `quality` includes end-of-run outcomes: use it for
eligibility/audit, never as a time-local predictive feature. Actual future control
events are also recorded for audit; select only events at or before the frame when
forming a feedback-policy input.

Evidence stays `EXPLORATORY`. Numerical PASS is separate from arrival/holding and
does not elevate physical evidence. Zero-flow and affine paired variants share one
train group; the original Y-CFD target family is test. Validation is empty because
there are only two independent initial/target families. Do not split frames or
paired controls to manufacture a third independent partition.

`physical_dt_s` counts actual time steps, not PBD iterations. A stopped failed
interval is separately preserved in `FAILURE_PROPOSAL.json`; its attempted step
count is additional to the last accepted saved prefix. The current fixed graph
stop localizes a violation to one0.5ms saved interval, not an exact failing substep.
Saved residence and arrival times have0.5ms sampling uncertainty at their edges.
CUDA event time spans run orchestration and host-idle gaps, so it is not exclusive
GPU active-kernel time. Throughput is valid complete trajectories per reserved
single-GPU hour, and must be quoted with bubble count and simulated horizon.
