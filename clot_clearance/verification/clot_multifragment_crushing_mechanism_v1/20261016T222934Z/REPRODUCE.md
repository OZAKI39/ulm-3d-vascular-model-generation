# Reproduce the frozen uniform / correlated-field comparison

Task `clot_multifragment_crushing_mechanism_v1`, accepted run `20261016T222934Z` (UTC). Research date in Europe/Zurich: 2026-10-17.

The server project is `/workspace/y_clot/project`; activate `/workspace/y_clot/env.sh`. Prepend this run's `build` directory to `LD_LIBRARY_PATH`. Use one MPI rank and one BLAS/OpenMP thread. Never overwrite this accepted run or an interrupted case.

## Frozen inputs

Read `SPECIMEN_AND_MATERIAL_FREEZE.json`, `INDENTER_GEOMETRY.json`, `INDENTER_LOADING_PROTOCOL.json`, `HETEROGENEITY_FREEZE.json` and `ANALYSIS_CRITERIA_FREEZE.json` before interpreting results. The original specimen files are copied byte-for-byte from the preceding mechanical bracketing run. Only one seed and one motion protocol are used. `HETEROGENEITY_FIELD.npz` contains the original white-noise grid, sampled field, original IDs and every local threshold.

The native indenter is a finite sphere using the upstream `fix indent`. Its force is added after the custom BPM force-density-to-N conversion. The Python indenter module is an independent oracle and work quadrature, not a force callback. The tool acts on all material nodes by current distance, including physical recontact after detachment. The prescribed trajectory is a kinematic impact surrogate, without an inertial hammer mass or a bubble/ultrasound model.

The isolated native change stores bulk damage onset by original bond ID; fixed-damage elastic energy and forces are unchanged. `sc`, interface parameters, reference weighted volumes and all original bonds are unchanged. Restart V3 persists the per-bond field; V2 inputs remain readable. `BUILD_AUDIT.json` records the upstream commit, source hashes, compiler/link commands and reused object hashes. Both the bond implementation and its generated class factory are recompiled. The old native library remains unchanged.

## Fresh reproduction

Copy this analysis directory to a fresh task directory and set its `TASK` and `RUN` constants in `common.py` to new values. Preserve the seed, physical constants, original input source and trajectory. `prepare.py` refuses an existing freeze; `run_case.py` refuses an existing history. This separation is deliberate: changing the accepted run in place is not a reproduction.

Run `prepare.py`, then `build.py` in the existing isolated server environment. Building reuses SHA-verified upstream objects listed in the preceding V3 build audit; a GitHub-only snapshot does not contain these objects. To rebuild on another host, reconstruct the locked upstream source and recorded native patches first, with the same compiler/package configuration. Do not point the new build at an arbitrary system LAMMPS.

With the new library first on `LD_LIBRARY_PATH`, run `checks.py`. It checks only the newly introduced current-contact force/energy/work path and original-ID onset storage. One combined sorting/restart fixture checks the newly persisted field. It is a constitutive consistency fixture, not a third physical impact case.

Run the uniform and heterogeneous physical cases in that order, each in a new output directory:

```bash
python run_case.py --name uniform
python run_case.py --name heterogeneous --heterogeneous
```

Use foreground supervisor wrappers for sustained jobs, following `manage.py` and the archived supervisor configs. Set `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`. The heterogeneous case must use exactly the same trajectory and numerical rules; do not select a new seed based on uniform results. No extra physical amplitude is permitted unless the task's explicit criterion is met (both cases have no softening).

## Evidence and reduction

Each case stores the full accepted scalar history, every accepted state's component statistics, a deduplicated original-ID component registry, sparse irreversible damage changes, native logs, event FP64 states and 50 µs visualization states. No unaccepted candidate becomes a result frame. The initial pre-cut is never counted as new fracture; contact pairs never become material connections.

Run `analyze.py`, `figures.py`, `report.py`, then `validate_outputs.py` against the new output root. `OPEN_RESULTS.html` links eight figure groups, the full final projection, CSV and review. The work integral follows the accepted piecewise-linear path, split at exact contact entry/exit; the ledger includes recoverable indenter energy and the existing particle-contact activation offset. Numerical residual is retained as residual. There is no numerical damping.

The original preflight variable-syntax failure is retained in `checks_attempt01_variable_syntax_failure`; `PRECHECK_IMPLEMENTATION_RECOVERY.json` records the equivalent LAMMPS expression correction. No physical loading or material parameter was changed by that correction.

After final reduction, update the WSL research document only with supported findings, preserve its before/after SHA and diff, and create a delivery manifest. Build objects and native binaries remain on the server; the delivery records their SHA and provenance instead of copying redundant solver builds. The original accepted trajectories, sparse events and key states are included in the local delivery.
