# F10S active code map

Baseline: complete F10R delivery `9f1db411836e9a67dff1f77a344e5457674e6a81`.

Bounded cleanup review: Git tracked status clean; F10R delivery 289 hashes verified. Module-reference search shows its helpers retain active execution or historical reproduction callers. No proven safe deletion; nothing deleted. All F1–F10R results, failures, snapshots, fields and response tables remain byte-identical.

- F10S orchestration: `scripts/straight_tube_exp1/f10s/{prepare,run,remote,analyze,diagnostics,figures}.py`.
- Immutable input protocol: `configs/straight_tube_exp1/f10s/protocol.json` and `results/straight_tube_exp1/f10s/PREREGISTRATION.json`.
- Source schedule: F10R `navigation_v1/runs/X_PAIR_DT1/INPUT.json`; single-source amplitude from its `ONE_HOT_X_SEQUENCE`.
- Frozen physics: `src/mevo_mb/fast_sim/{water_sequential,sequential_hard_sphere,water_b2,f10g_contact,nearfield}.py`; existing GPU fields and walls.
- Frozen acoustic/control interface: `acoustics/water_disc_proxy.py`, `control/scheduled_pzt.py`; F10N coarse/fine caches and APECSS table.
- Metrics: F10O `analyze.py`, F10P `analyze.py`; S duration adapter preserves every numerical threshold and checks the entire requested window.
- C/D only: analytic water-disc basis, original APECSS/F2, local frozen-shape resistance diagnostic; reference-material preflight. No hardware communication.
