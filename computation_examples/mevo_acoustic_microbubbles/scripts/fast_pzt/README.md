# Fast PZT aggregation F1

This is an explicit, isolated NVIDIA Warp backend. It does not replace any Stage 4/5 solver.
The delivered model is `fast_pzt_f1_v2_balanced_mobility`, with **EXPLORATORY** physics.

## Install and inputs

Run from `computation_examples/mevo_acoustic_microbubbles` on a CUDA-capable machine:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-fast-pzt.lock.txt
```

The tested isolated environment uses Python 3.12.3, Warp 1.18.0 and an RTX 4090
with driver 595.84. No system driver change or old-environment modification is required.
The lock file includes the small verification dependencies. `python.sh` sets PYTHONPATH
and single-thread BLAS; a full installation of the original acoustic solver is unnecessary.
The optional `fast-pzt` project extra is for an existing compatible project environment;
do not combine it with the old `warp-wall` extra's incompatible version pin.

Use these exact files under `results/fast_pzt/f1/inputs`:

| File | Existing server source |
| --- | --- |
| `response.h5` | `/workspace/mevo_acoustic_microbubbles/results/stage4/acoustic/response_extension/20261002T144654140341Z/adaptive_response.h5` |
| `basis.h5` | `/workspace/mevo_acoustic_microbubbles/results/stage3/acoustic/20261002T114813_234502Z_water_2p5MHz/acoustic/basis.h5` |

The small response file is archived in Git. The 33 MB basis is indexed by hash and
local/server path in `ARTIFACT_MANIFEST.json`; copy it from the existing project or
server archive, then verify the SHA-256 before running. Never substitute a different
frequency basis. The completed server deployment is `/workspace/mevo_fast_pzt_f1`,
with immutable executed sources under `.execution/<commit>`; `current` points to
the delivered source. Start new server runs from `/workspace/mevo_fast_pzt_f1/current`.

## Run a new dataset

Each label must be new. Existing output directories are never overwritten.

```bash
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt.experiment --label my_formation --config configs/fast_pzt/scene_b.json
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt.experiment --label my_off --config configs/fast_pzt/scene_b_off.json
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt.experiment --label my_pzt --config configs/fast_pzt/scene_c.json
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt.experiment --label my_batch --config configs/fast_pzt/batch8.json
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt.basic --label my_gpu_checks --checks-only
```

All delivered scene configs explicitly select 20 µs, one physical substep, and
12 local contact iterations. Keep the original numerical-quality criteria in
`f1_protocol.json`. Failure to pass them excludes a dataset from qualified generation.
Batch members are independent random clouds sharing radius, particle count and
the control schedule; there are no interactions between different batch members.

For server supervision, commit executable changes first, then use, for example:

```bash
python scripts/fast_pzt/remote.py experiment --label my_server_run --timeout 180 -- --config configs/fast_pzt/scene_b.json
```

This uses the existing SSH alias `vast4090` and creates a uniquely named supervisor
job. Scientific timing jobs should run one at a time. Check its `SUMMARY.json` and
logs; a process exiting alone is not a quality gate. Reproduce the final campaign with
`basic --checks-only`, `step_study --contact-iterations 12 --candidates-us 1 20`,
`performance --studies <step_label>`, `morphodynamics --performance <performance_label>`,
`throughput --performance <performance_label>`, and `final_scenes --performance <performance_label>`.
Every module takes a new `--label`. There is no need to rerun the retired prototype.

## Data contract and limits

Treat `trajectory.h5`, `SUMMARY.json`, `IDENTITY.json` and `CONFIG.json` as one dataset.
The HDF5 stores actual sampled positions, stable particle IDs, radii, physical time,
control indices, model/input identity, component memberships and shape tensors.
`SUMMARY.json` carries the completion and numerical-quality status, including
individual batch-member quality. A complete trajectory can still fail quality.
Each accepted-step adjacency change has a time and particle IDs in `edge_events`.
`TOPOLOGY.json` reconstructs simultaneous edge changes and verifies every saved
partition; its times are discrete accepted-step observations, not exact continuous
crossing times. Principal axes are unreliable near degenerate eigenvalues; use the
stored `axis_identifiable` flag. Cluster membership is a spatial graph, not adhesion.

Valid response domain: radius 1.779–1.781 µm, frequency 2.499–2.501 MHz,
pressure amplitude 0.1–40000 Pa. Delivered runs use 1.78 µm and 2.5 MHz.
The actual acoustic basis supports only 2.5 MHz. Zero drive has zero response;
nonzero pressure outside the table, frequency mismatch and out-of-grid positions
are errors, without clipping/extrapolation. No walls, background flow, fusion,
cluster locking, arbitrary velocity cap, global QP or coupled radial scattering
are included. Near-field regularization and local hydrodynamics remain uncalibrated.

## Figures and movie

On a plotting machine with Matplotlib, PyVista, Pillow and imageio-ffmpeg installed,
restore the indexed final HDF5 files and run:

```bash
PYTHONPATH=src:. python -m scripts.fast_pzt.analyze
PYTHONPATH=src:. python -m scripts.fast_pzt.visualize
```

These generate four English PNG/PDF figures with source data and a continuous
60 fps MP4. Coordinates and radii receive the same unit conversion; no cluster
rescaling, trajectory interpolation, synthetic motion or counterfactual stitching
is used. The camera is fixed. Each video frame represents 0.2 ms physical time.
`VIDEO_VALIDATION.json` verifies decoded frame count and frame rate.

Read `docs/FAST_PZT_F1_DEVELOPMENT.md` and `results/fast_pzt/f1/FINAL_GATE.json`
for measured performance and the boundary between numerical checks and physics validation.
