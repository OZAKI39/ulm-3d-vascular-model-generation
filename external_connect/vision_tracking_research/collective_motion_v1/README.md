# F10U-Vision / Collective Motion Perception

An executed, bounded observation prototype, built on F10T-Vision-2. **No persistent bubble IDs and no hardware actions.**
Development baseline `35755720a01de6ac6d6f14ac38e2a3ca84260236`; branch `dev/f10u-collective-motion-20261017`.

Start with [FINAL_RECOMMENDATION.md](FINAL_RECOMMENDATION.md) (Chinese 23-question report), [ACTIVE_CODE_INDEX.md](ACTIVE_CODE_INDEX.md), and [FINAL_GATE.json](FINAL_GATE.json).
Five PNG/PDF figures are in `figures/`; two fully decoded MP4 files in `videos/`. All estimates and failures in `results/`. The verbose real-observation JSON is tracked losslessly as `results/real_observation.json.gz`; a fresh checkout can restore it with `gzip -dk results/real_observation.json.gz`.

| Method | Endpoint tests | All-six output fraction | Accepted failures | Failure / accepted truth outputs | Non-affine correct rejection |
| --- | --- | --- | --- | --- | --- |
| PIV_multipass | 216 | 0.6991 | 21 | 0.1429 | 0.6667 |
| CPD | 216 | 0.5417 | 12 | 0.1081 | 0.5 |

This is not a general scientific accuracy pass. PIV nominal confidence intervals are undercovered; fixed texture/dropout cause false accepts. Real recordings have no six-parameter ground truth. `training_ready=false`, `LIVE_CAMERA_PASS=false`.

## Reproduce in the independent environment

```bash
cd /home/lzy/projects/external_connect/vision_tracking_research/collective_motion_v1
source .venv/bin/activate
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m unittest discover -s tests -v
python scripts/benchmark.py development
python scripts/benchmark.py test
python scripts/simulation_bridge.py
python scripts/real_observation.py offline
HF_HOME="$PWD/models/hf_cache" HF_HUB_DISABLE_XET=1 python scripts/ai_comparison.py
python scripts/window_response.py
python scripts/real_observation.py replay
python scripts/real_observation.py replay_fast
python scripts/supplement_checks.py
python scripts/analyze.py
python scripts/visualize.py plots
python scripts/visualize.py videos
python scripts/write_reports.py
```

Use a new output directory/copy for re-execution to preserve this delivered evidence. `real_observation.py camera` performs only a bounded existing Windows preview receiver test; it is deliberately not in the reproduction sequence above.
Dependencies: `requirements-lock.txt`. PyTorch/CUDA files were independently copied from the locally verified installation after a slow duplicate download; see `evidence/torch_dependency_copy.json`. No symlinks to the old environment and no changes to it. `pip check` passed. Official SEA-RAFT clone/weights are local and ignored by Git; see AI report for pins/downloads.
