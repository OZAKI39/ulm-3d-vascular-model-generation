# Fast PZT F2 reproduction

F2 extends the existing Warp engine. The selected exploratory backend is
`f2_near_field_candidate` + `f1_pbd_contact`, **10 µs**, one physical substep,
two local mobility iterations and twelve central PBD contact passes.
The archived nominal protocol retains 20 µs: its F2 overlap failure explains the
explicit 10 µs selection in delivered configs. The mobility-aware contact candidate
is available only for diagnostics; it failed the angular-balance audit.

## Environment and inputs

Use Python 3.12.3 with the unchanged `requirements-fast-pzt.lock.txt` (Warp 1.18.0).
The tested machine was an RTX 4090 with driver 595.84. From this project directory:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-fast-pzt.lock.txt
mkdir -p results/fast_pzt/f2/inputs
ln -s ../../f1/inputs/response.h5 results/fast_pzt/f2/inputs/response.h5
ln -s ../../f1/inputs/basis.h5 results/fast_pzt/f2/inputs/basis.h5
```

Only create absent links; reuse existing exact files. Restore the F1 basis from the
server path in `results/fast_pzt/f2/INPUT_MANIFEST.json` and verify SHA-256. The
small response file is already archived under F1. The 33 MB acoustic basis is not
copied into Git. Use only the verified 1.78 µm, 2.5 MHz configuration. The response
pressure domain is 0.1–40000 Pa; coupled effective pressures are also checked.

The isolated executed project is `/workspace/mevo_fast_pzt_f2` on `vast4090`.
It reads the existing F1 environment without changing it. Immutable source trees
are under `.execution/<execution-commit>`; `current` selects delivered source.
The final archive commit and historical execution commits have separate identities.

## Run new datasets

Every output label must be new; existing runs are never overwritten.

```bash
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.experiment --label my_formation --config configs/fast_pzt/f2/scene_b_adopted.json
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.experiment --label my_pzt_phase --config configs/fast_pzt/f2/mature_phase.json
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.experiment --label my_pzt_amplitude --config configs/fast_pzt/f2/mature_amplitude.json
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.experiment --label my_unchanged --config configs/fast_pzt/f2/mature_unchanged.json
```

The three mature continuation configs include the exact archived 32 positions and
an absolute-time offset of 102.2 ms; no random reinitialization occurs. The saved
mature state is complete for this overdamped memoryless, fixed-radius integrator;
forces/velocities are recomputed and no inertial/history state exists. Floating-point
atomic ordering means repeat trajectories are not bitwise deterministic.
`batch8_adopted.json` preserves the existing independent batching interface; the
F2 campaign tested only tiny batch independence, not a new eight-sample throughput.

Commit source changes before server deployment; do not run timing jobs concurrently:

```bash
python scripts/fast_pzt_f2/remote.py experiment --label my_server_formation --timeout 180 -- --config configs/fast_pzt/f2/scene_b_adopted.json
```

For the ordered F2 experiment sequence, use the following modules with new labels:

```bash
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.checks --label new_core
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.sensitivity --label new_sensitivity
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.adopt --label new_adoption --sensitivity new_sensitivity
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.mature --label new_mature --adoption new_adoption
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.phase_completion --label new_phase --mature new_mature
bash scripts/fast_pzt/python.sh -m scripts.fast_pzt_f2.performance --label new_performance --adoption new_adoption
```

Supervise server jobs with the wrapper and finite timeout. The underlying run has
a 120 s wall budget. Review `SUMMARY.json` quality, not merely process exit.
The initial sensitivity script's numerical selection is superseded by `adopt`,
which independently rejects the contact candidate's angular defect.

## Data and offline reproduction

`trajectory.h5`, `CONFIG.json`, `SUMMARY.json`, and `IDENTITY.json` form one dataset.
HDF5 contains all IDs, actual sampled positions, radii, control index, fixed field
identity, raw accepted-step pair-edge changes and cluster metrics. The quality
array retains cumulative maxima/corrections at each sample. Physical macro/substep
counts and sample overhead are explicit. `MORPHOLOGY.npz` separates translation,
Kabsch rotation and aligned internal deformation for fixed IDs and dynamic membership.
`TOPOLOGY_PERSISTENCE.json` never replaces raw edge events; 2 ms is the dwell rule,
not a physical bubble fusion criterion. `EDGE_DWELLS.npz` separately records brief
pair-threshold episodes, including censoring.

With Matplotlib, PyVista, Pillow and imageio-ffmpeg on a plotting machine, regenerate
the delivered `_v1` analysis and visualization after restoring indexed files:

```bash
PYTHONPATH=src:. python -m scripts.fast_pzt_f2.analyze
PYTHONPATH=src:. python -m scripts.fast_pzt_f2.validate_analysis
PYTHONPATH=src:. python -m scripts.fast_pzt_f2.visualize
PYTHONPATH=src:. python -m scripts.fast_pzt_f2.figures
PYTHONPATH=src:. python -m scripts.fast_pzt_f2.delivery
```

Analysis is explicitly wired to the delivered run labels; use a separate output
root or adapt the mappings for a new campaign. Do not overwrite original evidence.
Rendering changes no physics or particle positions. The MP4 contains 1112 observed
states at 60 fps (18.533 s), with uniform 0.2 ms physical spacing over 0–222.2 ms.
It follows the actual shared precursor and its exact-state amplitude continuation;
the equal restart state is not duplicated. All spheres have true 1.78 µm radius and
the camera/scale stay fixed. `STATE_PROVENANCE.json` and HDF5 source indices audit
each frame. No synthetic motion, cluster enlargement, interpolation or particle
subsampling is used.

Read `docs/FAST_PZT_F2_DEVELOPMENT.md` and the capability-specific `FINAL_GATE.json`.
The endpoint shape status is negative; transient eigenvalue responses are retained
separately. Contact sensitivity and lack of calibration must accompany reuse.
