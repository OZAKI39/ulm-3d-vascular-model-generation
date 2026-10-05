# Stage 2.5 reproducibility and guarded entry points

Use the existing server environment:

```bash
ssh -p 4159 root@50.115.148.16
cd /workspace/mevo_acoustic_microbubbles
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Branch: `stage2.5-pair-physics`. Read [STAGE2_5_VALIDATION_REPORT.md](STAGE2_5_VALIDATION_REPORT.md)
before interpreting outputs. Stage1 tag and Stage2 branch remain fixed.

## Completed evidence

Read `results/stage2_5/{zhang_forensic,response_table,pair_asymptotics,translation_map}/STATE.json`
and `results/stage2_5/audit/{PHYSICS_GATES,SOFTWARE_CHECKS,REGRESSION_STATE,CACHE_AUDIT}.json`.
The latest table is numerically verified under the original Tu/Marmottant assumptions,
not an experimentally calibrated universal SonoVue database:

```python
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
model = AdaptiveResponseTableBubbleModel(
    'results/stage2_5/response_table/20261001T220715485503Z/adaptive_response.h5'
)
response = model.query(1.78e-6, 2.5e6, 1.0e4)
# harmonic_response accepts complex pressure and applies h*phase per harmonic.
```

Fixed domain: R=1.6–1.96 µm, f=2.3–2.7 MHz, P=5–40 kPa.
Out-of-domain, nonfinite, incomplete and rejected-cell queries raise exceptions.
`backend='cupy'` provides bounded-memory batch queries.

## Gate check, with no simulation

```bash
.venv/bin/python -m mevo_mb gate-check run-lazarus --scope lazarus_paper
```

This currently exits nonzero with `NOT_RUN_GATE`. It is an admission check,
not a completed Lazarus simulator. The paper model, validated cylinder and
frequency sweep have not been unlocked. Nonlinear `bubbles run` checks gates
before reading the acoustic basis. `--force-demo` permanently records DEMO
metadata and watermark; it does not bypass numerical/response guards.
No override launched a simulation in this stage.

Classical transport can set `physics_validity.regime_map_path` to
`results/stage2_5/pair_asymptotics/STATE.json`. Warning thresholds are
`max_invalid_pair_fraction` (default0) and `max_transition_pair_fraction` (default0.05).
The map supports exact sampled uncoated uniform-drive points only. Unknown
parameters/shells or unequal drive phases remain outside.
`pair_backend: auto` uses the recorded crossover, or `pair_crossover` may be
specified explicitly. GPU unavailability/OOM records a NumPy fallback.

## Checks and managed reproduction

```bash
.venv/bin/python -m pytest -q
.venv/bin/ruff check src tests scripts benchmarks references
```

For future intentional long reproductions use supervisor, not a detached shell:

```bash
.venv/bin/python scripts/stage2_5/submit.py forensic_new scripts/stage2_5/zhang_forensic.py
.venv/bin/python scripts/stage2_5/submit.py asymptotics_new scripts/stage2_5/pair_asymptotics.py
.venv/bin/python scripts/stage2_5/submit.py adaptive_new scripts/stage2_5/adaptive_database.py -- --profile production --workers 7
```

Programs are named `mevo_stage25_<name>`, use `autostart=false`, do not automatically
restart, and log to `/var/log/portal`. Do not run simultaneous7-worker pools
under the7.68-core quota. These examples do not schedule work automatically.

`adaptive_database.py --resume <directory> --profile production --workers 7`
restores an identical build configuration. An extension creates a new child run:

```bash
.venv/bin/python scripts/stage2_5/submit.py extend_new scripts/stage2_5/extend_adaptive.py -- results/stage2_5/response_table/20261001T215645218867Z --max-nodes 24000 --max-depth 20 --seed 20261004
```

Resume an interrupted extension with the same parent/budgets/seed plus
`--resume <extension-directory>`. If failed holdout points inform further refinement,
use a new child and fresh holdout seed; those old points are no longer independent.

`zhang_forensic.py`, `pair_asymptotics.py`, `asymptotic_extension.py`,
`translation_map.py` and `benchmarks/contact_candidates.py` accept `--resume`.
Per-case checkpoints are reused. Reference resume refuses changed experiment
config or recorded physics-source hashes; do not bypass this check. Historical
development runs may differ from final source hashes, requiring a new run.

`ResponseNodeCache` shares deterministic keys across builds. Completed/failed keys
are not silently rerun; corrupt payloads fail loudly. Concurrent identical requests
are serialized. The first failed adaptive table retains `TABLE_REJECTION.json`;
never remove it to make an old query pass.

Execution CONFIG contains the then-current baseline commit and per-file source
hashes. Final source commit/result inventory are sealed separately. Raw APECSS
text is losslessly compressed on the server; local copies include node metadata
and retained cycles while raw `.txt.gz` remains server-side.
