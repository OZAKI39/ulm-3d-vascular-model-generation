# F10J execution map

This is a completed, bounded fixed-state investigation. Production integration is stopped because the admitted GPU implementation has no B=1 cost advantage.

Read `docs/FAST_PZT_F10J_DEVELOPMENT.md`, `FAST_PZT_F10J_GPU_SOLVER_AUDIT.md`, and `FAST_PZT_F10J_PERFORMANCE.md`. The layered decision is `results/fast_pzt/f10j/analysis/FINAL_GATE.json`.

Executed immutable remote labels, in order:

1. `audit_v1`: original unchanged reproduction and read-only PCG/Armijo trace.
2. `corrected_v1`: preregistered stable merit, 256 maximum PCG, original tolerances and explicit true residual admission.
3. `optimized_v2`: conditional captured PCG and accepted-state residual reuse; all 19 outputs match corrected positions and forces exactly.
4. `performance_v1`: one warmup and three paired measurements per implementation, plus one separate diagnostic profile.
5. `process_v1`: one new process per implementation, existing compiler cache.
6. `guards_v1`: expected incomplete-PCG rejection and actual GPU merit checks against extended precision.

These labels are immutable and should not be reused. `remote.py` only deploys committed source and refuses an occupied GPU. No source snapshot or historical result is edited on the server.

Offline reconstruction uses the existing shared WSL environment:

```bash
PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f10j.analyze_audit
PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f10j.topology_screen
PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.fast_pzt_f10j.figures
```

`archive.py` verifies each remote file and each source snapshot against its actual Git archive. `delivery.py`, `progress_document.py`, and `verify.py` complete reports, the minimal backed-up progress entry and final consistency checks. The progress edit refuses unexpected concurrent document changes. After any figure regeneration, visually inspect all four figures and update the explicit review hashes before running final verification.

No trajectory, B64 workload, new CFD, new acoustic basis, physical calibration or AI training is part of F10J.
