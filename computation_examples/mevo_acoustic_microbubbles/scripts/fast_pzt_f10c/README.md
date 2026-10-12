# F10C full-particle execution

Use the inherited Python/Warp environment. No F1–F10B module or dataset is replaced.
The production numerical candidate is `FullParticleBatch` with stride1 external fields.
`MultirateBatch` remains an **unadmitted research experiment** for ON-control morphology.

## Recorded campaign

- `profile_v1`: existing Graph comparison and bounded Warp activity timing.
- `batch_v1`: real serial/batch B1/4/16; `batch64_v1`: conditionally admitted B64.
- `multirate_v1`: pure decomposition and Q0/Q1/Q2/Q3 OFF/S/O.
- `checks_v1`: retained container-adapter assertion failure before physical checks.
- `checks_v2`: scope, injected failure attribution, cache/control/dt guards, memory audit.
- `failure_v1`: first Q2-S rejected proposal and bitwise accepted-prefix replay.
- `cold_v1`: actual fresh processes, identical outputs and sampled memory peaks.

The original batch/cold performance rows overlapped another project. They remain
shared-load observations. `exclusive_v1` waited for that project to exit, then replayed
serial/batch and fresh-process measurements with continuous sampled resource evidence.
Its B16/B64 results pass the 2x engineering gate; see the performance report for
measured values and scope. Do not interrupt unrelated jobs.

## New scientific job

Commit source/configuration first; select a new output label. The deploy tool refuses
to start if any CUDA compute process is already present (including other projects).
The existing Supervisor rules and finite timeout still apply.

```bash
python scripts/fast_pzt_f10c/remote.py profile --label profile_new --timeout 120
python scripts/fast_pzt_f10c/remote.py batch --label batch_new --timeout 240 --max-batch 16
```

B64 admission in the delivered script deliberately refers to the recorded B16 gate.
For a new campaign, register its gate/input identity before changing the reference.
`control_families` partitions records with independent sample IDs and explicit split
families by exact shared PZT schedule. It does not support asynchronous controls in
one batch. A device error aborts all members; do not label the others as valid.

## Local reconstruction and audit

```bash
PYTHONPATH=src:. python -m scripts.fast_pzt_f10c.interface_checks
PYTHONPATH=src:. python -m scripts.fast_pzt_f10c.analyze
PYTHONPATH=src:. python -m scripts.fast_pzt_f10c.figures
PYTHONPATH=src:. python -m scripts.fast_pzt_f10c.archive
PYTHONPATH=src:. python -m scripts.fast_pzt_f10c.verify
```

Use the shared WSL scientific Python for plotting. Figures use actual saved frames,
equal spatial scale, English labels, PNG and PDF. Original input hashes are reused;
no additional CFD, acoustic field or response table is generated. Source commit and
remote byte hashes accompany each dataset. All data remain EXPLORATORY and
`training_ready=false`; benchmark copies and numerical/control variants share one
training/test split family.

Exclusive replay (new label, idle device required):

```bash
python scripts/fast_pzt_f10c/remote.py exclusive --label exclusive_new --timeout 240
```

The driver stops its own child and rejects timing if sampled contention is found.
It runs at most one CUDA process at a time; the external monitor observes compute
PID cardinality and science Supervisor tasks, not pure GPU kernel duration.
