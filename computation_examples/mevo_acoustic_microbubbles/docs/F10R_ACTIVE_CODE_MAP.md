# F10R active entrypoints

Public baseline `467b17411deeb13a77370565b7382f071b41a774`; complete local science `41ae8ed2b88f2cbbed834fffccaa0b5de6e91589`. All 950 source/config/test/dependency files verified byte-identical. Different Git packaging, same scientific source. Historical root SYNC_MANIFEST.json is Stage 5; current packaging inventory is sync/20261014/.

No safe obsolete-code deletion was found in the bounded reference audit. No physics implementation or historical evidence was removed.

- New simulation-only protocol: src/mevo_mb/control/scheduled_pzt.py.
- F10R orchestration, diagnostics, checks and figures: scripts/straight_tube_exp1/f10r/.
- Frozen numerical reference: fast_sim/water_sequential.py, sequential_hard_sphere.py.
- Frozen fields/APECSS/F2/walls: water_disc_proxy.py, cubic_field.py, models.py, nearfield.py, water_b2.py, pdms_cavity_engine.py.
- Existing layered admission: scripts/straight_tube_exp1/f10p/analyze.py; fixed-ID metrics in f10o/analyze.py.
- Original inputs: results/straight_tube_exp1/f10n/.
- F10Q controls/trajectories: results/straight_tube_exp1/f10q/.
- New results: results/straight_tube_exp1/f10r/.

All times are ideal scheduled simulation events. No instrument transport or voltage conversion.

Final review entrypoints:

- `docs/FAST_PZT_F10R_DEVELOPMENT.md`: Chinese delivery with all 22 answers.
- `docs/F10R_SIMULATION_CONTROL_PROTOCOL.md`: units, clock, modes and immutable inputs.
- `docs/F10R_X_NAVIGATION_ROBUSTNESS.md`, `docs/F10R_Y_SHAPE_STABILITY.md`: actual measurements and limits.
- `scripts/straight_tube_exp1/f10r/analyze.py`: fixed-ID metrics and inherited layered admission; read control name from immutable INPUT.
- `scripts/straight_tube_exp1/f10r/verify.py`: independent all-step geometry, event pressure, clock, source and returned-file checks.
- `scripts/straight_tube_exp1/f10r/figures.py`, `report.py`, `update_progress.py`: data-only outputs and insertion-only Chinese history update.
- `results/straight_tube_exp1/f10r/analysis/FINAL_GATE.json`: PASS/FAIL/NOT_RUN layers, not a hardware approval.
- `results/straight_tube_exp1/f10r/analysis/SOURCE_METADATA_REVIEW.json`: first-run RESULT.source was overwritten by an output-loop variable after integration; raw results are preserved, analysis uses correct INPUT.source. Physics and controls were unchanged. New runner fixes the metadata variable.
- `results/straight_tube_exp1/f10r/diagnostics/AXIS_PAIR_PRESSURE_MAXIMA_REVIEWED.json`: continuous pressure-peak branches. Original Y nearest-center slope is explicitly invalid, retained in the original file.
- `results/straight_tube_exp1/f10r/DELIVERY_MANIFEST.json`: final file/hash inventory and source identity; `local/REMOTE_FINAL_STATE.json` confirms all four jobs exited and immutable execution source hashes remain unchanged.

Use the local shared environment with `PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1`; the remote environment remains `/workspace/mevo_acoustic_microbubbles/activate.sh`. Remote F10R root is `/workspace/mevo_acoustic_microbubbles/f10r_scheduled_control`; each run has a distinct `.execution/<commit>` and Supervisor label. Do not restart existing labels or overwrite preregistration. Results are complete; no further physics task is running.
