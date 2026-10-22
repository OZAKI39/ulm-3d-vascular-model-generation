# REALTIME_MICROBUBBLE_VISION_RESEARCH

Start with [中文选型结论](FINAL_RECOMMENDATION.md), [实测结果](REALTIME_TRACKING_BENCHMARK.md), [现有接口审计](EXISTING_INTERFACES_AUDIT.md), [开源项目调查](OPEN_SOURCE_TRACKING_SURVEY.md) and [未来接口设计](INTERFACE_DESIGN.md).

Live connection follow-up: [Chronos 2026-10-16](LIVE_CAMERA_FOLLOWUP_20261016.md). Actual Windows receiver verified at approximately 60 fps, including a separate run of the original launcher. WSL HTTP snapshots work; native WSL RTSP/UDP remains unresolved. `scripts/probe_live_camera.py` provides bounded read-only live checks, separate from the saved-video tracking benchmarks.

This is a bounded research harness, not a production tracking/control framework. **No hardware actions are implemented here.** The existing camera and AFG projects are protected and remain unchanged. Actual source video stays in `../sample_videos` and is not duplicated or uploaded.

## Environment

Independent `.venv`: Python 3.12; exact installed versions in `requirements.lock.txt`, package checks and GPU/OpenCV build details in `evidence/`. Official source clones and online checkpoint are in `third_party/`. See the survey for licenses; CoTracker has noncommercial terms and is not assumed BSD/MIT.

```bash
cd /home/lzy/projects/external_connect/vision_tracking_research
source .venv/bin/activate
export PYTHONDONTWRITEBYTECODE=1
python -m pip check
```

The older camera environment uses NumPy 2.x and Windows Python 3.14. Do not install this research dependency set into that environment.

## Reproduce without overwriting previous experiments

The original suite commands are recorded in logs. Each run refuses to overwrite an existing result directory. For a new replay, choose a new result tag and use the same library adapter:

```bash
python - <<'PY'
import sys
sys.path.insert(0, 'scripts')
from replay_benchmark import CFG, run
CFG['roi_xywh'] = [650, 350, 180, 350]
run('cc_norfair', 120, 210, 'independent_repeat_001')
PY
```

This reads only the saved video, inherits the Chronos latest-frame interface, and emits to a local mock consumer. It does not open camera or VISA hardware. Change the tag for each repeat; do not delete old runs merely to reuse a name. Input PTS cadence is unchanged. A large ROI is not equivalent to a small ROI performance result.

Derived results can be regenerated from existing canonical run files:

```bash
python scripts/analyze_results.py
python scripts/write_reports.py
python scripts/plot_research.py
python scripts/check_interfaces.py
python scripts/verify_delivery.py
```

These overwrite only derived research metrics/plots/check records, not raw source videos or original projects. Regenerating artifacts changes their hashes; regenerate the manifest when intentionally making a new delivery. Original execution commits are retained in local Git; this task does not publish a branch.

## Results layout

- `results/cc_norfair_*`, `blob_norfair_*`, `*_lk_*`, `*_farneback_*`: paired 90-frame clips, distinct broad/target ROI conditions.
- `results/full_cc_target_roi/`: first full 937-frame target replay.
- `results/quiet_confirmation_cc_target_roi/`: full post-installation confirmation, retained separately.
- `results/fault_cc_target_roi/`: slow consumer, simulated outage and timestamp discontinuity on real pixels.
- `results/cotracker_online/`: bounded 40-frame online-only GPU contrast and window ages.
- Each classical run saves `summary.json`, per-frame timings, decoder records, actual mock receipts, transient events, optional flow vectors and `states.jsonl`.
- `evidence/visual_annotations.json`: selected assistant-reviewed image landmarks; not expert biological ground truth.
- `evidence/interface_checks.json`: actual receiver reconnection and stale-message software checks, with no hardware.
- `figures/`: four actual-data PNG/PDF pairs. The architecture in `INTERFACE_DESIGN.md` is explicitly DESIGN.
- `INPUT_OUTPUT_SHA256.json`: scientific/research input/output manifest, excluding installed environments, downloaded third-party trees, `.git` and itself.

Finite diagnostic records are intentionally retained. They are not an unbounded live production logger. If later synchronizing this work, keep code, reports, figures, manifests and compact logs; do not add `.venv`, checkpoint weights or entire external repositories to the project Git branch.
