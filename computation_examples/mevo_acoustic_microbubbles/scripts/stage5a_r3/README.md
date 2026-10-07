# Stage 5A-R3 reproduction and review

Read `docs/STAGE5A_R3_PLAN.md`, `docs/STAGE5A_R3_ROOT_CAUSE_ANALYSIS.md`, and the
final machine-readable gate before interpreting motion as an accuracy result.
The project extends R2; no inputs or earlier result directories are regenerated.

## Execution identity

- Actual R2 archive base: `07db8ed35f3a460a4ba315e42b15d257f8518b0a`.
- Branch: `dev/stage5a-r3-wall-accuracy-performance`.
- V1 formal executable: `fa2fba5608037ed45143c8ede501cadf415e9180`.
- V2 formal executable: `b5079ad19be81ffe9e3d4f5b7c548717b3d43d67`.
- Server: `/workspace/mevo_acoustic_microbubbles_stage5a_r3`.
- Formal immutable trees: `.execution/<execution SHA>`.
- Outputs: `results/stage5/motion_solver_r3` (shared by immutable execution trees).
- Dependencies: `requirements-stage5a-r3.lock.txt`; environment and package-license
  inventories are in `results/stage5/motion_solver_r3/environment`.

The final Git archive SHA and the formal execution SHA have different meanings.
The validated R3 time-integration entry point is `scripts.stage5a_r3.run`.
Registering `translation_backend=constrained_lubrication_r3` supplies the R3 QP
adapter but does not by itself select this standalone R3 integration protocol in
an unrelated production driver. No production/reference promotion is implied.
Do not run `bootstrap.py` against developed R3 sources: it is the once-only
pre-development R2 identity check and intentionally asserts baseline hashes.
Do not overwrite an existing frozen protocol, holdout selection or campaign.

## Bounded local reproduction

On the existing server environment, from an isolated execution tree with the
verified R2 and Stage4L inputs available:

```bash
bash scripts/stage5a_r3/python.sh -m scripts.stage5a_r3.reproduce
bash scripts/stage5a_r3/python.sh -m scripts.stage5a_r3.wall_audit
bash scripts/stage5a_r3/python.sh -m scripts.stage5a_r3.run \
  --case B2 --duration 0.00001 --optimized --global-qp --pair cupy \
  --fixture results/stage5/motion_solver_r3/diagnosis/reproduction_v1/nominal/terminal.npz \
  --max-dt 0.000005 --tol-factor 0.25 --initial-dt 0.00000001 \
  --nonlinear-strategy quasi_newton --nonsmooth-recovery \
  --s0-dt 0.00000005 --budget 300 --label review_local_terminal
```

The local run is a full-system fixture continuation, not an original 3 s start.
The raw candidate NPZs contain full positions and QP arrays; `RESULT.json`
identifies exact finite features, raw gaps, continuous bounds and classifications.

## Source-bound tests and formal protocol

```bash
bash scripts/stage5a_r3/python.sh -m scripts.stage5a_r3.freeze_v2 --verify
bash scripts/stage5a_r3/python.sh -m pytest -q \
  tests/test_stage5a_r2.py tests/test_stage5a_r2_entry.py \
  tests/test_stage5a_r3_geometry.py tests/test_stage5a_r3_qp.py \
  tests/test_stage4l_lubrication.py tests/test_stage2_translation.py \
  tests/test_stage5a_r3_accounting.py tests/test_stage5a_r3_reporting.py
```

V2 nominal settings are B2 global QP / B3 component serial, one worker, BLAS=1, CuPy pair force,
quasi-Newton S1, bounded S0, initial dt 10 ns, maximum dt 5 µs, tolerance factor
0.25, S0 cap 50 ns, and output clock 10 µs. Refined settings are maximum dt
2.5 µs, tolerance factor 0.0625, and independently reduced S0 cap 25 ns.
All runs start from the original t0=3 s snapshots. Native calls and subprocesses
are watchdog-bound. The finite group and per-run budgets are sealed in the
protocol. A new reproduction campaign requires a new version/output directory;
the drivers refuse to overwrite their `CAMPAIGN_R3_V1` / `CAMPAIGN_R3_V2` directories.
V2 reuses the completed B2 V1 results with their original identities and carries
all prior charges forward. V1's untouched original scheduler ledger is resolved
by `CAMPAIGN_R3_V1/CLOSED_RESULT.json`; it is not an ongoing job. V2 follows a
separate source-bound 124-test gate and six newly selected numerical holdouts.

`PROTECTED_REPEAT_PLAN.json` separately preregisters two more protected B3 R2
+5 ms repeats within the existing 5400 s R2 group cap. The original V2 scheduler
ledger remains intact; `current.py` merges the separately closed supplement for
reporting and does not declare overall completion while it is pending. The
immutable follow-through executable is `248623f9ea7242d8ed938ef191217c98f036a625`.
Its managed watcher only reads completion state until all formal work closes,
then runs the supplement and six finite diagnostic windows sequentially.

## Postprocessing

Postprocessing can run locally after syncing the completed artifacts. Set
`PYTHONPATH=src:.` and use an environment with NumPy/SciPy/HDF5/Matplotlib/PyVista
and imageio-ffmpeg. `paths.py` translates preserved server references without
altering source identities or trajectory coordinates.

```bash
python -m scripts.stage5a_r3.resource_report
python -m scripts.stage5a_r3.formal_analysis
python -m scripts.stage5a_r3.reference_prefix --case B2
python -m scripts.stage5a_r3.reference_prefix --case B3
python -m scripts.stage5a_r3.analyze_early
python -m scripts.stage5a_r3.gate
python -m scripts.stage5a_r3.figures
python -m scripts.stage5a_r3.formal_figures
python -m scripts.stage5a_r3.early_figures
python -m scripts.stage5a_r3.event_video
python -m scripts.stage5a_r3.full_video --case B2
python -m scripts.stage5a_r3.full_video --case B3
python -m scripts.stage5a_r3.budget_ledger
python -m scripts.stage5a_r3.report
python -m scripts.stage5a_r3.review
```

Use `artifacts --server` only after scientific jobs finish, synchronize the server
manifest, then run `artifacts` locally. It validates strict JSON, inventories
large assets and original reference paths/hashes, and writes the size-bounded
`GIT_ARTIFACTS.txt` selection. The review directory retains full B2 fixtures and
complete lossless failed-run ledgers. `sync_completed` transfers only closed
run files with low-priority bounded I/O and records the transfer timing conditions.
After all managed scientific work closes, use `sync_completed --closed-delivery`
to freshly verify the campaign, protected repeats and diagnostic pipeline records
before allowing 4096 KiB/s transfers. The default remains 512 KiB/s for interim
copies. Neither mode reads growing run files or subtracts transfer overhead from
formal timing.

`dynamic_audit --window initial|stable|S0|event|B3_global_costly|B3_serial_costly` is an additional **post-campaign**
instrumented diagnostic. It refuses to run while formal timing is active. Its
exclusive cost partition, full/half1/half2 records and native-candidate records
are separate from completion timing. The accounting wrappers preserve solver
return values and exceptions; they do not modify scientific kernel files.
The extra B3 windows replay the actual +1.5 ms full-system state at the 5 µs
formal macro horizon to account for full/half1/half2 work at the costly state.
The review package retains complete losslessly compressed diagnostic JSON in
`review/compact_dynamics/`, including individual QP/native candidates and
full/half-step records. Decompress a file to its original `diagnosis/` name to
reproduce G7 directly; `DYNAMIC_COST_PROFILE.json` records the original hashes.

`early_precision` is the separately budgeted six-run diagnostic around the
first B2 resolution-target crossing. It uses each actual nominal/refined
+30 µs state and three levels, with explicit cold-cache initialization. It
does not replace an original-start formal run, and cross-lineage starts differ.
`B2_EARLY_PRECISION_AUDIT.json` retains this distinction and all 186 IDs.

Videos use physical sphere radii, identical old/new cameras, a uniform requested
physical clock and last actual stored states. They do not interpolate accepted
positions. Every encoded frame is decoded and its timestamp checked with the
available FFmpeg binary. Repeated actual-state frames are distinguished from
dropped or irregularly timed video frames. The local displacement plot has its
own explicitly labeled nanometre axes.

## Review priorities

1. `final/B2_EVENT_ROOT_CAUSE.json` and `final/WALL_CERTIFICATE_AUDIT.json`.
2. Full B2 nominal/refined fixtures and key accepted/rejected evidence in the
   lightweight review package; retain original source hashes.
3. `final/REFERENCE_LINEAGE_AUDIT.json`, all-active comparison NPZs, worst-ID
   neighborhood paths and same-clock sampling bounds.
4. `final/S0_S1_USAGE_AND_REFINEMENT.json` and `final/FORMAL_RUN_AUDIT.json`.

5. `final/DYNAMIC_COST_PROFILE.json`, `final/RESOURCE_LAYOUT_BENCHMARK.json` and
   full process wall times, excluding incomplete references from speedup.
6. `final/MOTION_SOLVER_R3_GATE.json`, source-bound JUnit, holdouts, manifest,
   source provenance correction and final delivery identity.

For +5 ms repetitions, all-active comparisons retain both repeat-0 versus later
R3 runs and each later R3 run versus its matching frozen/protected reference
repeat when available. Each gate row links only comparisons containing that exact
run; evidence is never copied from another repeat merely because the case matches.

No figure or smooth animation supersedes the quantitative precision gate.

Common-clock comparisons never interpolate positions. At unmatched legacy
times the neighboring-frame displacement is a scale, not a rigorous bound for
unrecorded projection motion. A certification flag now withholds a consistency
PASS when that guarantee is absent; original historical records are preserved
with an explicit clarification. Worst-neighbor figures use every stored
accepted state, with all physical axes on the same scale.
