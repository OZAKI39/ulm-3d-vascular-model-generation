# Asset identities and recovery

Critical source banks (~1.2 MB), exact scalar labels/index metadata, model-freeze identity, per-sample prediction arrays and reports are retained in Git. Large training image caches, raw recording, installed environments and neural weights are kept as local/server assets. No unselected third-party repositories are included.

Raw video SHA-256: `3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c`.

- Local raw recording: `/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4`.
- Server raw copy: `/workspace/f10v_vision_20261017/setup/f10v_source.mp4.complete` (verified identical).
- Local `training/<model>_seed<seed>/best.pt` and `last.pt`: exact hashes in each `summary.json` and final manifest.
- Server all checkpoints: `/workspace/f10v_vision_20261017/runs/prototype_7aed/training/`. Periodic epoch checkpoints are retained there in addition to best/last.
- Full cached data: server `/workspace/f10v_vision_20261017/data/`; local module `data/`. Regenerate using tracked source banks via `scripts/prepare_data.py --splits ...` without rereading raw video.
- Cross-environment tiny numerical differences are measured in `evidence/canonical_label_comparison.json`; original local label generation is retained under ignored `data/local_generation_labels/`, while canonical server labels are the formal authority.
- Exact server core archive: `/workspace/f10v_vision_20261017/setup/delivery_core.tar.gz`, archive SHA in `evidence/result_transfer.json`.
- Supplemental demo inputs and saved predictions: `results/synthetic_demo/`, with manifest hashes. Original F10U generator and fixed `configs/demo_protocol.json` allow recomputation only with frozen model weights.

Large JSON text records may be committed as `.json.gz`; run `python scripts/restore_text_assets.py` in a fresh checkout before plotting/analysis. This restores byte-identical originals and refuses to overwrite different existing content.

For an independent environment, install the pinned core requirements and the recorded official PyTorch/TorchVision CUDA build in a new venv. Consult `requirements-server.lock` and `evidence/dependencies.json` for all versions, source commits available from installed wheels and license files. Source commit information unavailable from a wheel is explicitly NOT_AVAILABLE, not invented.

No F10V weights are declared reliable or approved for physical control. Restore files for audit/research, keeping the frozen reject-all admission gate intact.
