# Stage2 working entry points

Local root: `/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles`.
Server root: `/workspace/mevo_acoustic_microbubbles`; SSH root@50.115.148.16 port4159 (`vast4090`).
Use server `.venv/bin/python`. Local base Python does not provide the scientific runtime.

Read [STAGE2_VALIDATION_REPORT.md](STAGE2_VALIDATION_REPORT.md) first. Stage2G has a retained failure;
H/J and the frequency sweep were not run. This page does not authorize resuming stopped BraVa/old mouse jobs.

```bash
cd /workspace/mevo_acoustic_microbubbles
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest tests -q
.venv/bin/ruff check src tests scripts benchmarks
```

Reproducible tasks, each creates a fresh timestamped directory:

| Command after `.venv/bin/python` | Action |
|---|---|
| `scripts/stage2/apecss_build_validate.py` | Build/test pinned official APECSS source |
| `scripts/stage2/radial_database.py --profile coarse --workers auto` | Broad table and holdouts |
| `scripts/stage2/radial_database.py --profile focused --workers auto` | Focused table; only passing cells available |
| `scripts/stage2/translation_validate.py` | Coupled N=1/force comparisons |
| `scripts/stage2/pair_validate.py` | N=2 references; preserves failure status |
| `benchmarks/pair_scaling.py` | Bounded synthetic performance, no aggregation |
| `scripts/stage2/kwave_benchmark.py` | Version-specific C++/CuPy fixtures |
| `scripts/stage2/regression.py` | Hash old results; recombine cached fields; rerun migration |

Set OMP_NUM_THREADS/OPENBLAS_NUM_THREADS/MKL_NUM_THREADS to1 for worker pools.
Long jobs on this server use the existing project-owned supervisor programs (`mevo_stage2_*`), not ad hoc detached shells.
They are autostart=false/autorestart=false; no final scientific job is left running.

If rebuilding elsewhere, clone APECSS from https://github.com/polycfd/apecss and checkout
`a0ee071ca247b5489e997d90e83b5f53c40e4ec7` under `third_party/src/apecss` before the build script.
Do not replace the pin with an unreviewed development branch. License/source records are in THIRD_PARTY.md.

`configs/stage2/strict_response_table.yaml` is an **unrun configuration template**, not a validated cylindrical case.
The focused table covers R1.6–1.96µm, f2.3–2.7MHz, p5–40kPa with holes in valid cells.
Initial radii, positions and local field pressure must be set inside that domain. Its water medium differs from old blood-like acoustic caches;
cache identity checks correctly reject cross-medium reuse. Physical translation viscosity and effective radial damping are separate inputs.
Production queries use `ResponseTableBubbleModel.evaluate_batch`/`harmonic_response`, never the external ODE inside a slow particle step.

`near_field_policy=strict` saves rejection reasons/state; `legacy_demo` labels NOT PHYSICALLY VALIDATED;
`pair_reference` permits only an explicitly selected uncoated uniform Gaussian N=2 problem with ZeroFlowField and no external acoustic field.
FlowField.velocity(x,t) is unchanged; generic VTKFlowField supports optional time snapshots with no space/time extrapolation.
