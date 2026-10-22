# Stage 5 handoff artifact index

## Recommended reading order for a new ChatGPT conversation

1. [docs/SYNC_STAGE5_20261006.md](docs/SYNC_STAGE5_20261006.md) — authoritative handoff and archive boundaries.
2. [Stage 4L final development report](docs/STAGE4_LUBRICATION_DEVELOPMENT.md) — completed work, incomplete runs and visualization history.
3. [docs/LUBRICATION_TRANSLATION_MODEL.md](docs/LUBRICATION_TRANSLATION_MODEL.md) — equations and approximation scope.
4. [Stage 4M FINAL_REPORT.md](results/stage4/secondary_ablation/final_decision/FINAL_REPORT.md) — failed A2 runtime gate, actual numbers and unanswered questions.
5. Current source: [translation backends](src/mevo_mb/bubbles/translation_backends.py), [tracker](src/mevo_mb/bubbles/tracker.py), [forces](src/mevo_mb/bubbles/forces.py), [secondary policies](src/mevo_mb/bubbles/secondary_ablation.py), [contact](src/mevo_mb/bubbles/contact.py), [cluster metrics](src/mevo_mb/bubbles/metrics.py), [VTK FlowField](src/mevo_mb/flow/fields.py).
6. [Stage 3 validation report](docs/STAGE3_VALIDATION_REPORT.md) — inherited validation scope, history-force and acoustic-grid limitations.
7. [Configs](configs/) and [tests](tests/), especially Stage 4L/4M fixtures and policy regressions.

## Source

- [src/](src/) — complete current scientific implementation, unchanged from the latest Stage 4M source state. Includes flow import, acoustics, radial responses, translation/contact and cluster infrastructure.
- [scripts/](scripts/) — complete current tooling, plus five previously untracked paper plotting scripts copied unchanged. Includes k-Wave basis tooling, APECSS workflows, Stage 4L/M runners and visualization tools. Archived scripts are not permission to run their campaigns.
- [pyproject.toml](pyproject.toml), [requirements.lock.txt](requirements.lock.txt), [third_party/](third_party/) — environment definitions and discovery/license metadata. External binaries, builds and environments are excluded.

## Tests

- [tests/](tests/) — complete test source and small dependencies; no test implementation was changed.
- [Stage 4L lubrication tests](tests/test_stage4l_lubrication.py), [Stage 4M policy tests](tests/test_stage4m_ablation.py) — 67 tests rerun during this sync; [JUnit evidence](results/stage4/sync_stage5_20261006/UNIT_TESTS.xml).
- [Sync verification](results/stage4/sync_stage5_20261006/VERIFICATION.json) — hash, scope, JSON/syntax/navigation and size checks; this is archive/software verification, not physical validation.

## Configs

- [configs/](configs/) — complete project configuration source. Historical run configs remain adjacent to each retained run identity.
- [Stage 4M policy documentation](docs/STAGE4M_NEAR_CONTACT_SECONDARY_ABLATION.md) — A0/A1/A2 definitions and frozen backend scope; full document.

## Stage 4L

- [Final gate](results/stage4/lubrication_upgrade/comparison_legacy/final_stage4l_20261005/GATE.json) — authoritative summary: `LUBRICATION_TRANSLATION_DEVELOPED_WITH_SCOPE`.
- [Secondary termination report](results/stage4/lubrication_upgrade/comparison_legacy/final_stage4l_20261005/SECONDARY_STAGE_REPORT.json) — incomplete budget-limited runs, preserved unchanged.
- [Original artifact index](results/stage4/lubrication_upgrade/comparison_legacy/final_stage4l_20261005/ARTIFACT_INDEX.json) — full original identities/hashes; some raw artifacts are intentionally outside Git.
- [Results collection](results/stage4/lubrication_upgrade/) — curated analytic checks, paired probes, primary baselines, timestep comparisons, diagnostics, figures, input snapshots and important failures. Summaries are complete copies; the raw collection is partial.

## Stage 4M

- [Final report](results/stage4/secondary_ablation/final_decision/FINAL_REPORT.md), [final decision](results/stage4/secondary_ablation/final_decision/FINAL_DECISION.json), [runtime gate](results/stage4/secondary_ablation/comparisons/M0_GATE.json) — complete original summaries; **A2_RUNTIME_GATE_FAILED**.
- [Dense comparison](results/stage4/secondary_ablation/comparisons/DENSE_COMPARISON.json) — A0/A1/A2 accepted steps, wall times, speedups and morphology diagnostics; complete summary.
- [Force-removal summary](results/stage4/secondary_ablation/comparisons/SECONDARY_REMOVAL_SUMMARY.json) — complete pair-norm ratios; these are not net force ratios.
- [L0 versus A0](results/stage4/secondary_ablation/l0_vs_a0/LUBRICATION_ONLY_FORMATION_EFFECT.json) — existing-data equal-time comparison at 0.08 s; not a new trajectory.
- [Dense run collection](results/stage4/secondary_ablation/dense_probe/) — exact configs/identities, A1 scope notes, initial snapshots, runtime/solver summaries and pre-integration failure records. Bulk per-step logs and trajectories are excluded.
- [Per-case diagnostics](results/stage4/secondary_ablation/comparisons/dense/) — small existing cluster time series and aggregate metrics; no new downsampling was generated.

## Historical validation

- [Stage 1 software report](docs/VALIDATION_REPORT.md), [Stage 2 report](docs/STAGE2_VALIDATION_REPORT.md), [Stage 2.5 report](docs/STAGE2_5_VALIDATION_REPORT.md), [Stage 3 report](docs/STAGE3_VALIDATION_REPORT.md) — full reports. Their scopes differ; Stage 1 is not a physical-validation certificate.
- [Stage 2 evidence](results/stage2/), [Stage 2.5 evidence](results/stage2_5/), [Stage 3 evidence](results/stage3/) — curated summaries, figures, configurations and small response-table references. APECSS cache nodes and full traces are excluded.
- [History-force sensitivity](docs/HISTORY_FORCE_SENSITIVITY.md), [known limitations](docs/KNOWN_LIMITATIONS.md), [physics provenance](docs/PHYSICS_PROVENANCE.md), [bibliography](REFERENCES.bib) — full interpretation and citation records; no uploaded copyrighted paper PDFs.
- [Sealed gold summary](results/stage4/performance_refactor/gold_baseline/STAGE4I_BASELINE_SUMMARY.md), [gold manifest](results/stage4/performance_refactor/gold_baseline/GOLD_BASELINE_MANIFEST.json), [backend decisions](results/stage4/performance_refactor/backend_decision/) — unchanged historical evidence; raw gold trajectories/logs are intentionally omitted.
- [Convergence](results/stage4/convergence/), [primary-only](results/stage4/primary_only/), [primary plus secondary](results/stage4/primary_secondary/) — curated identities, outcomes and summaries, preserving incomplete/failed states. Newer differing server metadata is captured separately under the sync evidence directory to avoid overwriting historical snapshots.

## Figures

- [M1: lubrication-only comparison](results/stage4/secondary_ablation/l0_vs_a0/Figure_M1_lubrication_only.png).
- [M2: dense performance](results/stage4/secondary_ablation/comparisons/Figure_M2_dense_performance.png).
- [M3: dense internal dynamics](results/stage4/secondary_ablation/comparisons/Figure_M3_dense_internal_dynamics.png).
- [M5: secondary removal](results/stage4/secondary_ablation/comparisons/Figure_M5_dense_secondary_removal.png).
- Corresponding PDF files are retained alongside these PNGs. M4/M6/M7 were not produced because the gates did not open.
- Stage 4L and historical validation figures remain beside their respective summaries. [Publication figures](results/stage4/publication_figures/) are existing exports; large duplicate SVGs and animation frames are excluded.
- [v11 interpretation note](docs/VISUALIZATION_V11_KNOWN_LIMITATIONS_20261006.md) — complete documentation of wall-sheet morphology and display-anchor flashes, without a rendering change.

## Reproducibility fixtures

- [tests/data/](tests/data/) — complete small unit fixtures, including four recorded contact-solver failures and the frozen Stage 4L force implementation.
- Stage 4L/M input snapshots, initial force/velocity arrays, source archives and run identities — exact copies, sufficient to identify short diagnostic inputs; full flow/acoustic files still require the omitted originals.
- Retained small adaptive response tables and response-database files are individually listed with hashes in the manifest. No new scientific data was generated for this archive.

## Previous archival provenance

- [sync/](sync/) — unchanged Oct-05 snapshot manifests, source-history bundle and project-specific orchestration logs. These are historical records; use `SYNC_MANIFEST.json` and the current handoff for this branch's contents and status. The original snapshot verifier targets the old export, not this curated tip.

## Excluded large raw artifacts

- [SYNC_MANIFEST.json](SYNC_MANIFEST.json) is the authoritative branch-content manifest: included-file sizes/hashes/origins, important results, environment and exclusions with server paths.
- [Excluded binary schemas](results/stage4/sync_stage5_20261006/EXCLUDED_HDF5_SCHEMAS.json) records HDF5 dataset metadata without loading trajectory arrays.
- Large HDF5/NPZ acoustic/trajectory arrays, huge event/solver logs, redundant source/graph/visualization dumps, caches and videos remain in their original local/server directories. Their absence is deliberate, not data loss or evidence of run failure.
- An original historical manifest can reference excluded files. Consult this index and the sync manifest before interpreting an absent dependency as missing source.
- Existing ancestor Git blobs are retained; there was no history rewrite, force push, main merge or original-result deletion.
