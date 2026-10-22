# Stage5A verification entry points

Status: **FAILED_WITH_EVIDENCE**. `constrained_lubrication_v1` is opt-in and has not passed the dense trajectory gate. Keep `lubrication_coupled` as the frozen reference. These files do not change acoustic/response/CFD or Stage4M A0/A1/A2 physics.

The earlier active-set campaign used `ce1808fd1` and was stopped by the user. A later explicitly authorized OSQP installation led to new, separate probe identities at `7dbc71f19f83687a8ff7e2646bc8402f1d5f634a`; both failed before accepting a step. The preregistered plan is from `b24fbc0fd`, with an OSQP supplement retaining all scientific tolerances. Per-run source/config/input hashes are authoritative. `SOURCE_COMMIT.txt` and `SOURCE_COMMIT_OSQP.txt` identify those executed revisions, not the final delivery commit.

## Inspect the delivery

- `docs/STAGE5A_MOTION_SOLVER_DEVELOPMENT.md`: engineering conclusions and scope.
- `results/stage5/motion_solver/OPEN_RESULTS.html`: eight PNG/PDF figures.
- `results/stage5/motion_solver/final/MOTION_SOLVER_GATE.json`: authoritative failed gate.
- `results/stage5/motion_solver/final/ATTEMPT_ANALYSIS.json`: rejection and solver-cost distributions.
- `results/stage5/motion_solver/diagnostics/`: B2 failed interval reconstruction and sparse matrices.
- `tests/data/stage5a_failed_qp/`: historical active-set failure fixture. It remains xfail for that backend and passes with OSQP. New OSQP dynamic failures have separate sparse fixtures in `diagnostics/`.

Raw baseline trajectories are read only from the server paths recorded in `SYNC_MANIFEST.json` and `audit_v2/BASELINE_SELECTION.json`; they are not duplicated here. The compact new accepted trajectories, failure logs, source hashes and plot arrays are included.

## Reproduce checks

Use project `/workspace/mevo_acoustic_microbubbles_stage5a`. The wrapper uses the existing base Python environment plus the installed `.solver_deps/osqp_1p1p3` overlay, with single-thread BLAS. Package versions and installation hashes are recorded; the base environment is unchanged. On another machine install `requirements-stage5a-osqp.lock.txt` into that overlay with `uv pip install --python BASE_PYTHON --target .solver_deps/osqp_1p1p3 --no-deps -r requirements-stage5a-osqp.lock.txt`, and provide the frozen project NumPy/SciPy environment. Set `MEVO_BASE_PYTHON` if its location differs.

```bash
bash scripts/stage5a/python_osqp.sh -m pytest -q tests/test_stage5a_osqp.py tests/test_stage5a_motion.py tests/test_stage5a_evidence.py tests/test_stage4l_lubrication.py tests/test_stage4l_metrics.py tests/test_stage4m_ablation.py
```

`compare`, `deliver`, `summarize_attempts` and `verify_delivery` are analysis-only modules. Figure rebuilding requires the recorded original baseline paths, the retained run summaries/HDF5 and the final JUnit XML. `verify_delivery` PASS concerns artifact integrity only; it never promotes the failed motion solver.

`audit_baseline --output NEW_DIRECTORY` measures independent original-state replays and refuses to overwrite an output directory. `benchmark` measures frozen states; its identical-state warm-start speedup is not a trajectory speedup.

The exclusive-run driver `run_probe` and sequential `campaign` enforce the joint 1 ms prerequisite for 5 ms. They are retained for reproducibility, but do not run them as an incidental verification step: the completed Stage5A campaign failed and no further scientific scope is authorized by this delivery. Future solver changes require a new plan and new run identities; never overwrite these failures or tune the frozen tolerances to obtain PASS.
