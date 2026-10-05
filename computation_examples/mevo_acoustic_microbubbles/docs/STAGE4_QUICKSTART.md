# Stage 4 execution and result inspection

Use the existing project at `/workspace/mevo_acoustic_microbubbles` on the server
(`ssh -p 4159 root@50.115.148.16`). The local mirror is
`/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles`.
The annotated `stage3-minimal-validated` tag and its accepted data are the baseline.

## Inspect the saved state

These commands are read-only and do not resume any calculation:

```bash
.venv/bin/python scripts/stage4/status.py
.venv/bin/python scripts/stage4/status.py --json
supervisorctl status mevo_stage4i5_candidates_20261005 mevo_stage4i5_completion_20261005
tail -n 25 /var/log/portal/mevo_stage4i5_candidates_20261005.log
tail -n 25 /var/log/portal/mevo_stage4i5_completion_20261005.log
```

The status script reads saved JSON and tolerates a worker progress file being
rewritten during the read. A recorded `RUNNING` status alone
does not prove that its process is alive; check supervisor and the age of the saved
progress. It deliberately avoids opening an HDF5 file while a worker writes it.
It also shows each child's `EXECUTION_FAILURE.json`, even while the parent waits
for another child in its process pool. A parent `RUNNING` record must not hide a
failed case. Incomplete cases remain excluded from scientific comparisons.
`results/stage4/CURRENT_EXECUTION_STATE.json` names the exact execution directory.
Each accepted component has its own `STATE.json`, pointing to an immutable run
directory with an input identity, source snapshot and summary.

## Protocol and provenance

`configs/stage4/bifurcation.yaml` defines geometry, steady flow, initialization and
seeds. `configs/stage4/experiment.yaml` defines integration and analysis. Formal
primary runs have 256 bubbles, ten paired seeds, a 6 s observation window and
0.5 ms base timestep. Residence time begins at the initialization plane x = −2 mm.
This is a finite initial ensemble, not the separate BraVa perfusion protocol.

Within each seed B0–B3 share the exact initial positions, radii and particle IDs.
Their identity hashes, CFD input, acoustic field and response table are recorded
in `RUN_IDENTITY.json`. Parent stages preserve `source_snapshot.zip`.
The controls change both phase and amplitude. A conditional exit fraction
`N_A / (N_A + N_B)` must be read alongside `N_A / N_initial`, `N_B / N_initial` and
the unassigned fraction. An unassigned bubble is censored at the observation
cutoff; it is not a successful outlet arrival or a numerical failure.

Stage 4 uses the existing SimVascular solver and native results. The optional
`tetra_p1` import path verifies tetrahedral membership before interpolation.
Production still excludes history forces, streaming and feedback. The validation
reference in `validation/history_branch_robustness` is an adiabatic two-scale
approximation with separately reported interpolation, contact and wall limits.

## Continue an interrupted ordered execution

Check the saved execution, actual process status, prerequisite identities and
completed case files first. Do not submit a second copy of a running stage.
Long calculations on this instance use the finite supervisor helper; jobs have
`autostart=false` and `autorestart=false`.

The original B2/B3 have naturally completed and are sealed under
`results/stage4/performance_refactor/gold_baseline`. **Never rerun the full gold.**
The old October 3 ordered supervisor is finished; it must not be restarted.
`run_remaining.py` now rejects H/I when that gold manifest exists.

The active continuation is documented in [STAGE4I5_PROTOCOL.md](STAGE4I5_PROTOCOL.md):
the candidate service finishes its matched short replays; the completion service
then runs the additional ordered-contact candidate, the correctness-first backend
decision, revised J, K, final audit and report. Both are finite, non-restarting
jobs. Do not submit another copy of either service. J/K refuse to start without
the completed I.5 backend decision. For a failed downstream phase, inspect and
reuse completed evidence before a new invocation; do not repeat earlier stages.

The original per-particle wall certificate and global NNLS stay available.
The compiler and Warp backends are optional, with recorded environment/binary
identities and independent regression gates. History stays outside production.

Choose a new descriptive job name for each invocation. Starting a stage creates a
new timestamped directory. It is not a command to resume an arbitrary partial
particle trajectory: accepted substep state in an incomplete HDF5 file is failure
evidence, not a complete restart checkpoint. Preserve failed attempts and only
reuse closed cases after checking the full identity and protocol. Never restore
the separately stopped old BraVa or mouse tasks as part of this workflow.

## Review and regenerate reports

After all required stages exist, the final audit runs regression, archived
reference tests, lint, Stage 3 file-integrity checks and Stage 4 input/result
identity checks. Finalization consumes those results and produces:

- `STAGE4_ACCEPTANCE.json` and `docs/STAGE4_VALIDATION_REPORT.md`;
- eight scientific figures, six performance figures and `OPEN_RESULTS.html` in the timestamped report run;
- native CFD, labelled geometry and active-particle VTP/PVD exports for ParaView;
- per-particle branch-switch, crossing and censoring tables, plus secondary and
  runtime summaries.

These paths describe generated deliverables; their existence and component
statuses determine whether the execution has finished. `WARNING` is a scientific
outcome, not permission to change a label to `PASS`. Numerical agreement does not
validate omitted near-wall hydrodynamics or establish experimental accuracy.

To regenerate a report from completed evidence, use a new
supervisor invocation starting at `audit` or `report` as appropriate. Finalization
creates a fresh report directory and updates the latest report/acceptance pointers;
it does not rerun particle trajectories or overwrite previous run directories.
