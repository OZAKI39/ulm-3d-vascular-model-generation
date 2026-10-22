# F10E bounded execution

All scientific runs were executed once in the registered order. Existing labels are immutable; do not rerun them for analysis or silently replace a failed result.

1. `isolation --label isolation_v1`: four original singleton retries, only after F10D identity audit.
2. `SELECTION.json`: earliest independent S/O failure plus original preparation1104.
3. `diagnostics --label contact_v1 --phase baseline`: two passive-window C0 replays.
4. `diagnostics --label alternatives_v1 --phase alternatives`: C1/C2/C3 to 5 ms only.
5. `EXTENSION_DECISION.json`: select the first passing affordable candidate under the registered rule.
6. `diagnostics --label extension_v1 --phase extension`: one original-state 10 ms replay, C1.

The remote wrapper deploys committed code to `.execution/<commit>`, checks the whole GPU, and creates a unique Supervisor program. It refuses existing labels or an occupied GPU. `timed` records full child-process wall time and sampled GPU PIDs. The first isolation job predates this high-resolution envelope; its Supervisor timestamp interval is retained explicitly.

Use `python -m scripts.fast_pzt_f10e.remote diagnostics --label <new-authorized-label> --phase <phase>` only when a new scientific run is actually warranted. A label change alone is not scientific justification for repeating these cases.

Offline commands, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src:.`:

```bash
python -m scripts.fast_pzt_f10e.analyze
python -m scripts.fast_pzt_f10e.figures
python -m scripts.fast_pzt_f10e.verify
```

`archive` audits fixed remote snapshots against recorded input/source hashes, transfers result bytes, and verifies the source against local Git objects. It never changes a scientific result or restarts a process.

`ContactAuditBatch` uses the F10B observer via multiple inheritance with F10C; no original kernel is replaced. It checks the original error after each physical step, records only the registered contact window, and exports the last 32 audited steps. No incomplete checkpoint is used. The runner's new step/class/observer arguments are opt-in; existing defaults are unchanged.

Each numerical variant has a separate mode ID but retains the same seed split group. Complete controls, preparation trajectories, short-window tests, retries, and diagnostic copies remain distinct categories. All are `EXPLORATORY / PHYSICS_UNCALIBRATED`, `training_ready=false`.
