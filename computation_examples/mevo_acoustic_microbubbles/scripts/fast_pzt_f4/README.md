# F4 reproduction

Start from branch `dev/fast-pzt-f4-y-geometry-20261008`, under `computation_examples/mevo_acoustic_microbubbles`. F1/F2/F3 source and results remain independently usable. Read `docs/FAST_PZT_F4_Y_GEOMETRY_DEVELOPMENT.md` and `results/fast_pzt/f4/FINAL_GATE.json` before interpreting the evidence.

## Inspect archived evidence locally

Use the existing shared Python environment (NumPy/SciPy/h5py/Matplotlib/PyVista/VTK/Pillow/imageio-ffmpeg). These commands require no Warp or GPU:

```bash
export PYTHONPATH="$PWD/src:$PWD"
python -m scripts.fast_pzt_f4.analyze
python -m scripts.fast_pzt_f4.visualize
python -m scripts.fast_pzt_f4.visualize --video
python -m scripts.fast_pzt_f4.validate
```

`read()` supports full HDF5 and the lossless `SAVED_STATE/METRICS/DIAGNOSTICS.npz` archive fallback. Every saved frame and raw edge event is retained, not temporally decimated. The figure/video source hashes identify the original HDF5 even when viewing from compact copies. Rendering units are converted consistently to mm globally and µm locally; relative bubble/wall sizes do not change.

## Prepare the exact original inputs (only for a new experiment checkout)

Restore `inputs/labelled_lumen_m.vtp`, `PLANES.json` and the reused `basis.h5`/`response.h5` from `INPUT_MANIFEST.json` and verify SHA-256. The archived prepared mesh and finite-radius initial states may be reused directly. To reproduce preparation, run both commands in order in an isolated checkout:

```bash
python -m scripts.fast_pzt_f4.prepare
python -m scripts.fast_pzt_f4.prepare_controls
```

The second command reproduces the documented pretrajectory field-based control/duration amendment; do not run only the first and mistake its initial plan for the final protocol. Commit prepared configuration before launching. Set `PYTHONPATH` as above.

## Run on the existing server

The existing server root is `/workspace/mevo_fast_pzt_f4`; its environment is linked to the F1 `.venv` with Warp 1.18.0. Do not overwrite existing labels. The launcher deploys committed source into `.execution/<SHA>` and creates a bounded, non-restarting supervisor job. Run stages sequentially and inspect their RESULT before the next; avoid GPU concurrency during performance measurement.

```bash
python -m scripts.fast_pzt_f4.remote checks --label checks_new --timeout 180
python -m scripts.fast_pzt_f4.remote aggregation --label aggregation_new --timeout 360
python -m scripts.fast_pzt_f4.remote controls --label controls_new --timeout 240 -- --aggregation aggregation_new
python -m scripts.fast_pzt_f4.remote edge_checks --label edge_checks_new --timeout 120
```

`performance.py` deliberately gates on the delivered `aggregation_v1` and `controls_v1` evidence and reuses F3's N200 input; it may be run with a new output label after reviewing those gates:

```bash
python -m scripts.fast_pzt_f4.remote performance --label performance_new --timeout 240
```

No command above resumes an old batch or Stage 5A-R3. Science uses real PZT fields, 600 ms formation and 120 ms continuation at 10 µs. Prescribed velocities appear only in `checks` / `edge_checks`, labelled GEOMETRY_UNIT_TEST. Performance N200 is uniform-pressure PERFORMANCE_ONLY. `vessel.enabled=false` uses original F3 semantics. `wall_contact_enabled=false` is a separate F4 switch that retains open-boundary handling while disabling physical WALL projection.

## Inputs and outputs

Authoritative geometry identity, field-overlap definition, validation, aggregation/control summaries and performance are separate root JSON files under `results/fast_pzt/f4`. `EXECUTION_INDEX.json` maps completed and failed attempts to source SHAs. `ARTIFACT_MANIFEST.json` excludes itself, identifies the archive policy and records local/server paths. `PRESERVATION.json` checks the inherited F1/F2/F3 manifests without altering them.

The final seven acceptance categories are deliberately separate. Numerical geometry PASS does not mean validated wall hydrodynamics or robust shape control. Physics scope remains EXPLORATORY / GEOMETRIC_WALL_ONLY.
