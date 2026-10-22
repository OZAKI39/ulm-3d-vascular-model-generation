# F10AC-QD + F10AD-Multifrequency

This extension reuses the F10AA CPU contact integrator and Warp acoustic kernels. It does not change `src/mevo_mb`, the historical 2.5 MHz `Segment`, or historical result files. Actions are Pa-equivalent, not voltage. All results remain `EXPLORATORY`, `PHYSICS_UNCALIBRATED`, and `training_ready=false`.

## Read the completed study

The authoritative entry is `results/straight_tube_exp1/f10ac_qd/FINAL_CHINESE_REPORT.md`. Read `FINAL_GATE.json` for scoped admission, `ELITE_MANIFEST.json` for the actual trajectory behind each archive cell, and `ALL_RUNS.json` / `AUDIT.json` for completed and rejected attempts. `COMPARISON_V2.json` is the corrected fair comparison. The first comparison is retained as an integration audit and included in cost; do not mix the two rankings.

The report's total includes measured scientific and postprocessing stages. It excludes development, package download/build, and transfer time. Charged comparison costs allocate the common eight evaluations equally to all algorithms; these are not three separate physical executions.

## Environment and execution identity

Server root: `/workspace/mevo_acoustic_microbubbles/f10ac_ad_20261018`.

Use its **own** `.venv`, not the shared earlier environment. Exact package versions, build script and third-party licenses are under `results/straight_tube_exp1/f10ac_qd/evidence/`. Python 3.12, Warp 1.18.0+cu12, SciPy 1.15.3, pyribs 0.12.0 and pinned APECSS are used. APECSS is an external executable under `third_party/build/ultrasound_apecss`.

Every executed source tree is an immutable Git archive in `execution/<commit>_<label>`. Its `SOURCE_F10AC_COMMIT.txt`, `run.sh`, Supervisor configuration and source archive hash are preserved. `evidence/server_logs/*_EXECUTION.json` gives each real command and commit. Report commits do not replace execution identities.

Within a verified snapshot, the common runtime environment is:

```bash
export PYTHONPATH=src:. MPLBACKEND=Agg
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv/bin/python -m pytest -q tests/f10aa tests/f10ab tests/f10ac
```

Historical inputs resolve through read-only project result symlinks. `common.verify()` checks the frozen configuration, initial library, historical caches/response table, and original source SHA-256 before scientific execution. Local `.venv` is only a convenience link to the previous local environment; it is not the authoritative dependency set.

## Executed stages and code

| Stage | Entry | Output |
|---|---|---|
| Freeze | `f10ac_qd.prepare` | Protocol, initial IDs, protected hashes |
| Regression and resources | `f10ac_qd.pilot` | Two bitwise historical trajectory regressions, phase gauge, 1/2/4 workers |
| Fresh fields and radial response | `f10ad_multifrequency.build` | 2.4/2.5/2.6 MHz FNM caches, 324 APECSS solves, separate field/response gates |
| Fair search | `f10ac_qd.search` | Official BOP-Elites, CMA-MAE, Sobol; raw suggestions and fixed archives |
| Bounded extensions and validation | `f10ac_qd.workflow` | Winner extension, discrete-frequency search, independent seeds, dt/cache checks |
| Retained out-of-bounds review | `f10ac_qd.outside_validate` | One preregistered candidate, three fresh seeds; failed convergence retained |
| Independent audit | `f10ac_qd.audit` | All accepted steps of complete trajectories; compact datasets and raw hashes |
| Interpretation/media | `f10ac_qd.interpret`, `figures`, `videos` | Conservative mode interpretation and actual saved-frame presentation |
| Delivery | `f10ac_qd.report`, `finalize` | Reports, per-frequency gates, old-result protection and SHA-256 manifest |

Commands above are Python module suffixes under `scripts.straight_tube_exp1`. Existing outputs are frozen research evidence. A new campaign requires a new output root and preregistration; do not rerun `prepare/build/search` over this completed study. `remote.deploy()` creates a new immutable source snapshot and a unique, bounded Supervisor job; its exact executed invocations are in the execution records.

## Data semantics

- Full HDF5: every accepted position step, actual attempted-step forces/gates and fixed particle IDs. Failed prefixes end at their actual physical time; no padding.
- Portable `dataset/*.h5`: exact saved frames every 10 steps plus endpoint, raw frame indices and original raw SHA, every attempted metric/gate, sampled forces and contact reactions. Self-contained, without external HDF5 links. These are not every-step archives and cannot replace all-step validation inputs.
- Local portable copies cover all video inputs and the three response-domain failures; a subset of full raw files is also present. All 366 full raw trajectories and all portable datasets remain on the server, indexed by SHA. Remote raw storage is part of the delivery. `LOCAL_AVAILABILITY.json` distinguishes local and remote files.
- For the three acoustic-query failures, the legacy raw `failure_proposal_m` may be the previous accepted proposal. The derived attribute and audit explicitly distinguish it from an actual rejected integrator step. The reported low pressure is an independent CPU reconstruction at the last accepted state.
- `METRICS.npz` and summaries retain 3D and XY shape statistics, Kabsch/affine diagnostics, distances, connectivity and force/contact quality. Geometric proximity is not physical adhesion or bubble fusion.
- The video reader accepts either full raw HDF5 or its portable HDF5 with verified parent hash and raw-frame indices. Videos show exact stored positions with fixed camera, equilibrium sphere radius and spacing. No physical-time interpolation or invented shape smoothing. 2 ms is displayed at 30 fps over 6.7 s; frame-to-raw mappings and decode checks are saved.

## Scope

The main three candidates pass their 24 numerical comparisons and three fresh-seed checks. The additional strong-contraction candidate fails three new-seed dt comparisons and is not an admitted new mode. The nine fixed-action frequency trajectories pass 18 dt/cache checks, within their listed initial state/action/2 ms scope. This does not admit all controls at those frequencies.

Two non-collinear endpoint descriptors across mixed frequencies/actions do not establish a rank-two local control Jacobian. One repeated sign category is confirmed; no resolved absolute principal extension was found. Y responses are `HOMOGENEOUS_WATER_PROXY_SPECIFIC`. No new experimental calibration, real PZT execution, arbitrary shape control, or world-model training is claimed.
