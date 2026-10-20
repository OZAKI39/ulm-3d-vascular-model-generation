# V4 reproduction

Use the isolated server environment `/workspace/y_clot/env.sh`, single MPI rank and one BLAS/OpenMP thread. Native library SHA: `4510aed9195a30b03a4895435be05eca005d462507a12c4423fe8d2a854481e2`. The original V3 baseline archive and V2 native build objects are required: this compact Git snapshot is not a standalone binary distribution.

1. Preserve existing accepted archives. Run `python code/restore_export.py` in a compact Git checkout to restore and verify the losslessly compressed native reference, body data and correspondence tables. Reproduction must use a fresh project copy or a new RUN configured in `code/common.py`; output guards deliberately refuse replacing a frozen reference or trajectory.
2. Verify `SOURCE_FREEZE.json`, `V4_IMPLEMENTATION_FREEZE.json`, environment package records and original V3 sources.
3. Run `run_p0.sh`, `run_p1.sh`, then `performance_summary.py` to reproduce the paired eight-step profiling, before generating the irregular reference. Do not repeat the historical full V3 dynamics.
4. Run `prepare_reference.py`, `build_native.py`, and `build_surface.py`. The native builder verifies and reuses the frozen V2 object hashes, recompiles changed classes/contact/dispatch, and writes isolated V4 outputs. Full commands and source/object hashes are in `BUILD_AUDIT.json`.
5. Prepend the new output `build` directory to `LD_LIBRARY_PATH` **before** importing LAMMPS. Run `check_elastic.py`, `check_native.py`, then `admit.py`. The tests use separate small validation coupons/fixtures; they are not additional formal cases.
6. `run_formal.sh` runs `run_case.py pilot` to 2.0ms and `run_case.py continue` from its legal checkpoint. There is one physical case. If a numerical gate fails, stop: do not weaken gates, refit materials or add expensive cases.
7. After terminal status, run `terminal_audit.py`, `analyze_result.py`, `render_videos.py`, `render_target_closeup.py`, `write_report.py`, `validate_delivery.py`. `reference_figures.py` and `compact_performance.py` prepare the reference-only outputs and small performance payload. `finish_visual_review.py` repeats only the revised spectrum plot, visual rendering, report and delivery checks, without dynamics or repeated scientific analysis.

Current service: `supervisorctl status clot_v4_formal_20261019T221000Z`. Log: `/workspace/y_clot/project/verification/clot_irregular_lps_discretization_v4/20261019T221000Z/formal.supervisor.log`. Postprocessing is a separate managed service, never another physical trajectory.

The initial generation uses exactly seed20261019, bounded isotropic ball proposals, fixed exterior/notch/anchor positions, a documented geometric strong-volume projection, actual mapped-cell volumes, and KD-tree original neighborhoods. Never substitute nominal equal volumes in a V4 replay.

Full original FP64 states, native restarts, build objects and dense damage streams remain at `/workspace/y_clot/project/verification/clot_irregular_lps_discretization_v4/20261019T221000Z` on the server. The Git delivery includes code, manifests, small scientific tables, plots and selected evidence; large source payloads are identified explicitly in export manifests. Videos and dense data must be copied from the result archive if excluded by the local Git size policy.
