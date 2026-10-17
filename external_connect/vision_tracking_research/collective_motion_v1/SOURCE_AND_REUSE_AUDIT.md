# Source and reuse audit

Vision repository: `/home/lzy/projects/external_connect/vision_tracking_research`. Baseline `35755720a01de6ac6d6f14ac38e2a3ca84260236`. Tracked baseline clean; numerous old untracked GUI evidence files preserved (`evidence/git_status_before.txt`). No applicable AGENTS.md in the workspace ancestors. No subagents used.

Read prior supervised report, frozen human masks/manifests, held-out 49/50 report, OpenCV and Norfair interfaces, state freshness/mock code, Windows bridge and original LatestFrameStream. OpenCV held-out center F1 0.862 applies only to previously reviewed ROIs, not full-frame boundaries or arbitrary experiments. No labels or detector thresholds were tuned here. No Cellpose/SAM/SAM2 training.

`src/distribution.py` calls the original `scripts/replay_benchmark.py::Detector('cc')`, threshold 10. It scopes/restores the in-memory ROI configuration. Norfair is a dependency of that legacy module, **no Tracker is instantiated**. Candidate quality is a contrast heuristic, not probability. Counts are visible candidates, not verified bubble identities.

Windows project `/mnt/f/CHRONOS_control`, WSL original `/home/lzy/projects/external_connect/CHRONOS_control`, and Tektronix project were located/read as applicable. Original Windows receiving code and `smart_online_tracking_v2/src/live_bridge.py` reused unchanged. No Tektronix module is imported by the new observation pipeline. Existing receiver made 3 built-in attempts within a single bounded invocation, received 0 frames. No network configuration changed.

True latest applicable simulator worktree: `/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10s_20261015/computation_examples/mevo_acoustic_microbubbles`; HEAD `31a44d3f342e5d0bca55f7f9b83f4d32b8189696`, branch `dev/fast-pzt-f10s-x-pair-onehot-trap-scope-20261015`, tracked clean. Separate public sync HEAD `c7363522ef492d09f7902c0c3382572a70beb68a` was not assumed to be the working research baseline. The older unsuffixed mevo directory is not this F10S worktree.

Input video `/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4`, SHA-256 `3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c`. 1476×1120, recorded PTS cadence 0.0333 s (~30.030 Hz); no microscope px/um or exposure timestamp calibration. Protected input hashes before/after are in `evidence/`; source and output hashes in `INPUT_OUTPUT_SHA256.json`.

Safe cleanup: **zero deletions**. No active temporary code was established to be unreferenced and without scientific value. No physics runs, hardware controls, LED tests, SAM2 or world-model training. New independent Python 3.12 environment only. Failure logs and the first OpenPIV masked-array API failure remain saved. The API fault was repaired using official masked arrays before independent testing.
