# Active code index

The preserved engine is `src/mevo_mb/fast_sim/accelerated_water.py` + `accelerated_sequential.py`, original response/force and cavity code. New work lives in `scripts/straight_tube_exp1/f10af/`:

| Script | Purpose |
|---|---|
| prepare.py, field.py, static_extra.py | Freeze identities; extended caches; force equilibria/winding/nulls |
| states.py, metrics.py | Real legal states; F10Y-compatible descriptors |
| run.py | Original F10AA Campaign.run with new inputs |
| slow.py | Restricted SciPy RK45 with original fallback |
| pilot.py, paired.py, extension.py | Replay, paired checks and time ladder |
| generalize.py, longtime.py, sobol.py | Independent tests; bounded extension; full-physics fallback |
| tooling.py, tooling_smac_fix.py | Official API replay, explicitly not formal physics search |
| audit.py, analyze.py | Accepted-step audit and scientific classification |
| figures.py, videos.py, report.py | Actual-data media and evidence-indexed reports |
| remote.py, final_analysis.py | Immutable snapshots, Supervisor, history/tests |
| rotation_checks.py, repair_failed_view.py | Exact shared-time motion checks; derived failed-prefix correction |
| receive.py, deliver.py, repair_environment.py | Read-only receipt; document/hash seal; isolated environment launchers |

Configuration: `configs/straight_tube_exp1/f10af/protocol.json`. Tests: `tests/f10af/test_metrics.py` plus inherited AA/AB/AE regressions. No duplicate particle, QP, GP or Hyperband solver was introduced.

## Evidence

- `PREREGISTRATION.json`
- `TEST_RESULTS.json`

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false
