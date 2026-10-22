# F10V-Vision — Geometry-Aware Collective Microbubble Motion

Actual trained prototype; **not scientifically admitted to replace PIV**. Read [Chinese findings](FINAL_RECOMMENDATION.md) and [development report](F10V_DEVELOPMENT_REPORT.md). No physical device actions. `training_ready=false` concerns future physical world-model data, not whether these vision models were trained.

Code baseline:9cb22d729fece92fe160082edf09222a18068d74. Training snapshot7aed42436; frozen evaluationed36aa804; continuous demoa5a4ef53d. Full identities and source hashes accompany artifacts. Branch:dev/f10v-geometry-aware-motion-20261017. Local commit only; no F10V GitHub push.

Local module: `/home/lzy/projects/computation_examples/_archives/mevo_vision_f10v_20261017/external_connect/vision_tracking_research/collective_motion_v2`. Remote: `/workspace/f10v_vision_20261017`, env `.venv`, data `data/`, output `runs/prototype_7aed/`, immutable code `snapshots/`. These folders are independent of the physics and old vision environments.

Quick review: `figures/Fig_F10V_1...7.*`, `videos/Video_F10V_1_real_causal_replay.mp4`, `Video_F10V_2_synthetic_truth.mp4`, `Video_F10V_3_difficult_ablation.mp4`. All plots are actual data; rejected raw diagnostics are clearly labeled.

Run in this module after activating the isolated environment:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/figures.py
.venv/bin/python scripts/videos.py
```

Full reproduction should use a **new output and cache directory** so frozen results are preserved. Restore the real-patch bank to the new cache; it is only about1.2MB and has per-file hashes. The raw video is needed only to re-extract banks or replay real microscopy. Regeneration without raw video:

```bash
export F10V_DATA_ROOT=/absolute/new/cache
export F10V_OUTPUT_ROOT=/absolute/new/run
# Copy data/banks and bank_manifest.json into F10V_DATA_ROOT first.
.venv/bin/python scripts/prepare_data.py --splits train validation calibration final_test --workers 8
.venv/bin/python scripts/train.py --all --batch 192
.venv/bin/python scripts/evaluate.py freeze
.venv/bin/python scripts/evaluate.py test
.venv/bin/python scripts/baselines.py --split final_test --workers 4
.venv/bin/python scripts/analyze.py
.venv/bin/python scripts/replay_runtime.py --video '/absolute/path/source.mp4'
```

The training script has a bounded budget and `--resume`; use the server Supervisor configuration retained in `evidence/server_execution/` for long runs. Never rerun `test` to tune a model; the script refuses to overwrite a completed test. Reproduction in a fresh output directory is an audit, not a new independent test. Canonical server labels and source hashes are retained; floating differences across environments are explicitly measured.

Offline inference:

```bash
.venv/bin/python scripts/infer.py --images frame0.png frame1.png frame2.png --pts 0 0.0333 0.0666 --model geo_fused --output /tmp/f10v_state.json
```

Current frozen thresholds reject all AI outputs. This is intentional evidence-based behavior; do not bypass it for hardware control. `--legacy` exports an explicit1.0-compatible message. All input images must be same ROI and use actual increasing timestamps.

Weights are retained locally under `training/*/best.pt` and `last.pt`; periodic checkpoints remain on the server. Weight hashes are in run summaries and `INPUT_OUTPUT_SHA256.json`. Virtual environments, full caches, third-party downloads and weights are excluded from Git; `ASSET_RESTORE.md` documents recovery. Some large text records are stored as `.json.gz` in Git and automatically restorable with the provided script.
