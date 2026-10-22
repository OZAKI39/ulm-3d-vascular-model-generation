# F3 reproduction and evidence

Run from `computation_examples/mevo_acoustic_microbubbles` on branch
`dev/fast-pzt-aggregation-f3-20261008`. Read `docs/FAST_PZT_F3_DEVELOPMENT.md`
and the frozen `configs/fast_pzt/f3/protocol.json` before interpreting results.

The server uses `/workspace/mevo_fast_pzt_f3`, its immutable `.execution/<SHA>`
source folders, and the existing `/workspace/mevo_fast_pzt_f1/.venv` (Python 3.12.3,
Warp 1.18.0, NumPy 2.2.6, SciPy 1.15.3, h5py 3.16.0; RTX 4090). The `current`
link denotes the final delivered source; it is not the execution identity of older
jobs. Each trajectory's `IDENTITY.json` is authoritative for its execution commit.
No server environment or original F1/F2 data is modified.

Inputs are the exact archived F2 mature state, F2 basis, and response table.
`results/fast_pzt/f3/INPUT_MANIFEST.json` supplies hashes and original locations.
On the existing server, F3 `inputs/basis.h5` and `inputs/response.h5` are read-only-use
symlinks to the preserved F2 inputs. A new checkout can fetch the basis from the
manifest server path; do not rebuild a replacement basis or response table.
F2 historical model contrasts must remain available for sensitivity analysis.

## Inspect or regenerate the archived analysis

No GPU is needed for the following commands. The local shared environment contains
NumPy, SciPy, h5py, Matplotlib, PyVista/VTK, Pillow and imageio-ffmpeg; the exact
versions are in `VISUALIZATION_ENVIRONMENT.json`.

```bash
export PYTHONPATH=src:.
export OPENBLAS_NUM_THREADS=1
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f3.analyze
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f3.validate
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f3.visualize
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f3.figures
```

These regenerate derived files from the closed `driving_v3`, `checks_v2`,
`repeats_v1`, `controls_v1`, and `performance_v1` results. They do not advance or
alter raw trajectories. The validator permits absent excluded historical F1 files
locally and reports them; separate server preservation validation covers all files.
The visualizer uses actual radii/positions and identical camera bounds. It precedes
`figures` because the latter embeds its actual 3D force panels.

## Re-run bounded GPU experiments

Use a unique run label. `remote.py` deploys committed code and starts a bounded,
non-restarting supervisor job on SSH alias `vast4090`. Never reuse an existing
label or delete old data to make a command run. For example, a new matched phase
trajectory using the full archived mature coordinates is:

```bash
python -m scripts.fast_pzt_f3.remote experiment --label phase_repro_v1 --timeout 180 -- --config configs/fast_pzt/f3/mature_phase.json
ssh vast4090 'supervisorctl status mevo_fast_pzt_f3_phase_repro_v1'
```

Other ready configurations are `mature_amplitude.json`, `mature_unchanged_stable.json`,
`mature_unchanged_atomic.json`, and `performance_N200_{atomic,stable}.json`.
Default F3 execution is `stable_reduction`; `fast_atomic` preserves F2 atomic sums.
`stable_pairs`/`stable_components` exist only for diagnostic ablations. All enforce
the same frozen physics parameters and preserve independent microbubble motion.

To reproduce the bounded scientific sequence, launch each line only after the
preceding job reports completion and its `RESULT.json` exists:

```bash
python -m scripts.fast_pzt_f3.remote driving --label driving_repro_v1 --timeout 180
python -m scripts.fast_pzt_f3.remote checks --label checks_repro_v1 --timeout 180
python -m scripts.fast_pzt_f3.remote repeatability --label repeats_repro_v1 --timeout 240 -- --checks checks_repro_v1 --driving driving_repro_v1
python -m scripts.fast_pzt_f3.remote controls --label controls_repro_v1 --timeout 240 -- --repeats repeats_repro_v1 --driving driving_repro_v1
python -m scripts.fast_pzt_f3.remote performance --label performance_repro_v1 --timeout 120 -- --repeats repeats_repro_v1
```

This is a reproduction recipe, not an instruction to resume computation on opening
the project. Each new job writes its own raw evidence/analysis. The top-level
archived figure/summary scripts intentionally continue to reference the original
closed F3 labels, so reproductions cannot silently replace the published comparison.

`delivery.py` assembles module-specific gates, input/source indices and the artifact
manifest after closed evidence has been synced. It does not claim a universal PASS
and does not run science. `FINAL_GATE.json` preserves transient shape response,
negative robust endpoint findings, field limitations and morphology sensitivity.
