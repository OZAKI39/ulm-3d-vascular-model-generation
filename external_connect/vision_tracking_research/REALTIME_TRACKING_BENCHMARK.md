# Real-time tracking benchmark — 2026-10-16

**CAUSAL_REPLAY_TEST. LIVE_CAMERA_NOT_AVAILABLE at the original probe. No physical controller was used.**

Subsequent 2026-10-16 read-only checks verified actual Windows camera reception at approximately 60 fps and WSL HTTP snapshots. Native WSL RTSP remains unresolved. See [live connection follow-up](LIVE_CAMERA_FOLLOWUP_20261016.md). All tracking performance below remains causal saved-video replay data, not live camera tracking.

## Main finding

The selected OpenCV + Norfair path is computationally practical in a limited ROI, but the reviewed identity stability is not yet adequate for automatic bubble control. Positive landmark detection and per-frame message throughput must not be mistaken for verified bubble identity or successful group tracking.

## Input and protocol

- Input: `../sample_videos/Recording 2026-10-15 185901 截取视频.mp4`, SHA-256 `3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c`. H.264, 1476×1120, 937 decoded frames. Consecutive PTS interval 0.0333 s; first/last PTS 0/31.1688 s; container duration 31.232 s. The container average-rate field differs slightly from PTS cadence. Use actual PTS, not an assumed 60 fps.
- Recording provenance beyond the user-provided sample is not independently verified. The filename suggests a recording/crop. Camera recording FPS, exposure time and scaling are UNKNOWN. No physical µm scale or physical µm/s is reported.
- Three common windows per method: frames [120,210), [360,450), [720,810), each 90 frames. Broad ROI [640,80,620,960], with fixed stripe x830–990 masked; target ROI [650,350,180,350], at original pixel scale. ROI conditions are kept separate.
- Detection calibration uses first frames 0/30; threshold/area/association parameters are frozen in `protocol.json`. The broad calibration failure is retained (`smoke_cc`, 662 initial unmasked components in calibration image). Smaller ROI was added after the broad smoke failed timing; it does not represent full-field capacity.
- OpenCV local background correction and connected components / SimpleBlobDetector; Norfair 2.3 Euclidean gate 8 px, initialization delay 2, hit counter 5. No bespoke tracking algorithm or learned detector training. Missing frame steps explicitly advance the filter; motion units use actual PTS.
- PyAV producer decodes sequentially at the original PTS cadence; no future images are supplied to the tracker and no backward trajectory smoothing is used. Initial file seek/preroll is excluded from steady-state timing. Original Chronos latest-slot consumer methods are reused.
- Latest decoded-frame slot is bounded to one. Consumer counts overwritten frames. Diagnostic result lists grow only within the bounded experiment and do not feed future observations back into processing. Output uses local ZeroMQ HWM=1/conflated receiver and sequence/age checks. The primary age limit is 100 ms.
- Quality-screened visual state = at least one fresh, observed, sufficiently persistent/moving/high-contrast/non-neighbor-ambiguous image candidate. It is **not** a verified biological target. All messages have hardware permission false and admitted bubble-control count zero.
- WSL2, Intel i7-13700HX (24 logical CPUs), Python3.12 isolated venv, OpenCV4.11.0.86/NumPy1.26.4/Norfair2.3.0, OpenCV one thread. GPU only used for the separate CoTracker test: RTX4060 Laptop 8 GiB, torch2.11.0+cu128.
- Shared laptop, not an exclusive laboratory timing rig. Dependency downloads/installation overlapped early comparisons; no simultaneous tracking benchmarks ran. A post-install full-replay confirmation is reported separately. CPU percent uses one full CPU core as 100%; process RSS includes retained finite diagnostic records. Thermal/OS/background activity contributes to tails.

## Paired method comparison

Each row pools actual measured per-frame values over the same three windows within that ROI. Quantiles are from raw frame timings, not averaged per-clip p95s. Rate denominator is measured wall time including startup/close; source PTS cadence is 30.030 Hz.

| ROI | Method | Processed / 270 | Processing mean / p95 ms | Drop % | Processing output Hz | Mock quality-screened visual Hz | CPU % | Max RSS MB |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| wide | cc_norfair | 178 | 50.55 / 70.20 | 34.07 | 18.94 | 16.49 | 123.2 | 267.2 |
| wide | blob_norfair | 159 | 57.79 / 91.18 | 41.11 | 16.81 | 11.95 | 118.2 | 276.9 |
| wide | cc_norfair_lk | 162 | 55.41 / 71.72 | 40.00 | 17.36 | 14.68 | 119.6 | 274.7 |
| wide | cc_norfair_farneback | 88 | 105.12 / 122.14 | 67.41 | 9.20 | 0.00 | 115.1 | 289.1 |
| target_roi | cc_norfair | 270 | 7.31 / 9.29 | 0.00 | 29.40 | 27.99 | 47.7 | 209.2 |
| target_roi | blob_norfair | 270 | 5.09 / 6.64 | 0.00 | 29.50 | 27.53 | 37.9 | 213.1 |
| target_roi | cc_norfair_lk | 270 | 9.59 / 11.81 | 0.00 | 29.43 | 28.01 | 52.5 | 220.9 |
| target_roi | cc_norfair_farneback | 270 | 16.69 / 20.98 | 0.00 | 29.38 | 27.96 | 74.7 | 225.0 |

The broad ROI failed to follow input cadence. Dense Farnebäck there yielded **zero fresh valid visual messages** under the admission limit; it must not be called a 30/60 Hz controller. Target-ROI Farnebäck followed this 30 Hz sample but p95 exceeded the 16.7 ms reference budget for 60 fps. Sparse LK is the lighter auxiliary. Blob was faster in the target ROI but missed more reviewed landmarks.

## Full primary replay and confirmation

| Run | Processed/source | Processing mean/p50/p95/max ms | Decoded / processed / valid visual wall Hz | Skipped |
|---|---|---|---|---|
| Full target-ROI replay | 935/937 | 7.67/7.82/9.51/39.11 | 29.97/29.91/29.72 | 2 |
| Post-install confirmation | 936/937 | 7.58/7.50/9.84/64.82 | 29.98/29.94/29.75 | 1 |

Confirmation was run because background dependency installation was an unresolved timing confound; it is not parameter tuning or replacement of the first result. Both raw runs remain available. Figures retain the first full run. No 60 fps camera/input test occurred.

### Per-stage measurements: full primary target ROI

| Stage | Mean ms | p50 ms | p95 ms | Max ms |
|---|---:|---:|---:|---:|
| decode_ms | 4.182 | 4.297 | 7.458 | 13.429 |
| frame_age_at_start_ms | 0.373 | 0.379 | 0.672 | 6.491 |
| detect_ms | 1.624 | 1.596 | 1.938 | 4.780 |
| track_ms | 2.519 | 2.570 | 3.186 | 6.649 |
| shape_quality_ms | 3.137 | 3.162 | 4.253 | 33.578 |
| flow_ms | 0.001 | 0.001 | 0.001 | 0.014 |
| serialize_send_ms | 0.388 | 0.390 | 0.524 | 1.144 |
| processing_ms | 7.669 | 7.820 | 9.509 | 39.115 |
| receive_to_send_ms | 8.042 | 8.115 | 9.953 | 39.815 |
| transfer_ms | 1.900 | 1.404 | 1.943 | 108.282 |
| receive_to_mock_ms | 9.557 | 9.218 | 11.264 | 114.049 |

`processing_ms` starts at consumer pickup and includes detection, tracking, per-region shape/quality, optional flow, serialization/send. `receive_to_mock_ms` is **SOFTWARE_PROCESSING_LATENCY from post-decode receipt**, including wait/transport/receiver deserialization. Decode time is separate. These overlapping intervals cannot be summed as independent stages. `CAMERA_TO_CONTROLLER_END_TO_END_LATENCY = NOT_MEASURED` because exposure/synchronization/live path are absent.

Broad pilot code timed skipped-filter advancement in the detection column; the saved `evidence/pilot_benchmark_source.py` preserves this exact implementation. Before target tests, accounting moved that work into tracking. Total processing times are comparable; do not compare the broad-versus-target detect-only subcolumn as if accounting were identical. Early pilot results predate the first execution commit and are identified by their archived source SHA; later runs record Git HEAD and source hash.

## Accuracy: what was and was not established

Manual inspection here means **assistant visual review**, not a trained experimentalist annotation or biological ground truth. `evidence/visual_annotations.json` records 21 selected clear positive landmarks in three keyframes and three short trajectories, ±3 px localization uncertainty. Annotations were made before reviewing comparison tracker outputs; they are positive-only, not exhaustive labels.

| Detector | Matches / selected positive landmarks | Misses in this selected set | Biological precision / recall |
|---|---:|---:|---|
| cc | 20/21 | 1 | NOT_EVALUATED |
| blob | 9/21 | 12 | NOT_EVALUATED |

Matching is one-to-one within 6 px. This selected-landmark hit rate is not full bubble recall. Unmatched detections may be real unannotated features, debris or artifacts; a false-positive rate cannot be derived from this label set.

| Online method | Landmark | Observed matches / reviewed timestamps | ID changes in matched observations | Mean position error px | Reviewed span s |
|---|---|---:|---:|---:|---:|
| cc_norfair | landmark_A | 1/4 | 0 | 0.50 | 0.0000 |
| cc_norfair | landmark_B | 2/4 | 1 | 0.88 | 0.1998 |
| cc_norfair | landmark_C | 1/3 | 0 | 0.53 | 0.0000 |
| blob_norfair | landmark_A | 1/4 | 0 | 0.43 | 0.0000 |
| blob_norfair | landmark_B | 1/4 | 0 | 1.52 | 0.0000 |
| blob_norfair | landmark_C | 1/3 | 0 | 0.44 | 0.0000 |
| cc_norfair_lk | landmark_A | 1/4 | 0 | 0.50 | 0.0000 |
| cc_norfair_lk | landmark_B | 2/4 | 1 | 0.88 | 0.1998 |
| cc_norfair_lk | landmark_C | 1/3 | 0 | 0.53 | 0.0000 |

First timestamps are absent while tracks initialize. Later reviewed landmarks also disappear or receive different IDs; the small position error **conditional on a match** cannot compensate for missing/unstable identity. The detected ID change is fragmentation evidence, not proof of a two-object identity swap. Reliable same-ID pairs were insufficient in the chosen short labels for a meaningful velocity-error estimate; **velocity accuracy remains NOT_EVALUATED**, not zero.

The first full replay produced 709 temporary software IDs and 621 LOST→OBSERVED transitions; these are diagnostic software events, not 621 confirmed physical re-identifications. Motion changes and image artifacts can generate false connections. Long-lived software track duration is not automatically legal identity duration.

Group morphology: the adapter actually computes current threshold-region area, contour, centroid, second-moment tensor, axis ratio and orientation. However, the sample has no reliable labeled aggregate segmentation. Consequently `group_shape=null`, `group_shape_status=UNKNOWN_UNVALIDATED_SEGMENTATION`; do not report these component shapes as validated microbubble group shape. Small/large aggregates, physical fusion, expert shape error, global false joins/ID switches, target-crossing error and true fluid velocity remain **NOT_EVALUATED**.

The real clip contains dense fine texture, motion, long fixed structures and changing defocus; close features are visible. A clean sparse-bubble case, verified physical exit/re-entry identities and verified local rotation are not independently annotated. No synthetic scenario is substituted for those missing cases.

## Fault injection and actual receiver checks

On real frames [120,240), an added 120 ms consumer delay, 350 ms producer outage and PTS jump were injected. These are **synthetic software faults on real images**, not an observed camera failure. 34/120 frames were processed; 86 overwritten; valid visual messages 0. The latest slot remained bounded to 1. Two timestamp discontinuities reset tracking. Received-frame age p95 159.00 ms; all stale states were rejected.

`scripts/check_interfaces.py` separately executes the original Chronos receiver loop with a hardware-free mock backend returning actual saved image pixels: receive failure clears the slot, reconnect resumes frames, stop exits cleanly. A real ZeroMQ mock test verifies stale, duplicate, discontinuous and untrusted messages are rejected. These are software tests, **not live reconnect or network-loss validation**. Results: `evidence/interface_checks.json`.

## CoTracker3 Online bounded GPU contrast

Official source `82e02e8029753ad4ef13cf06be7f4fc5facdda4d`, official scaled-online checkpoint hash in `evidence/cotracker_weights.json`. Input: same target ROI and real frames [120,160), 9 initial image landmarks, internal resize [384, 512]. Maximum 16 queued image inputs, window16/step8, frame arrival cadence unchanged. No offline tracker call, no whole-video preprocessing.

The first 16-frame window spans 499.5 ms; nominal next step is 266.4 ms. Actual first-window inference was **802.23 ms** and the second inference **206.79 ms**. Newest frame ages at release were **824.57 ms** and **255.32 ms**. Both failed the 100 ms age gate; 0 admitted fresh emissions. The first overrun caused a new online session/feature seed rather than silently carrying IDs across a skipped window. This measures two bounded emissions only: the second is not a statistically established steady-state latency. All 9 points being visible does not establish correct bubble identity.

Peak CUDA allocated 1589.0 MB; reserved 1906.3 MB; GPU NVIDIA GeForce RTX 4060 Laptop GPU. No out-of-memory failure. Performance is specific to this local laptop and small experiment. It is not a fair full-pipeline accuracy comparison against a detector+tracker with automatic births; no claim of CoTracker superiority is made.

## NOT_RUN / failure ledger

- Bonsai: source/documentation/platform feasibility only; no Bonsai/.NET runtime installed in WSL, no Windows workflow executed. Local ZeroMQ test is not Bonsai validation.
- DeepLabCut-Live / Bonsai.Sleap: no task model; no heavy training or unrelated-model FPS claim.
- btrack: source/API audit only; no incremental-stream performance claim.
- GStreamer appsink: documentation/version audit only; existing FFmpeg path retained.
- Live camera: both bounded TCP probes timed out; previous 60 Hz session logs remain historical evidence only.
- Initial broad detection found fixed texture/blur fragments; broad pipelines dropped frames; target landmarks suffered identity instability; neural windows missed age deadline. These failures are retained.
- Some initially requested official source/doc paths returned 404 (e.g. Bonsai VideoCapture.cs); repository tree inspection located actual CameraCapture/FileCapture source. No missing endpoint was represented as a working API.

## Reproducibility and artifacts

All figures are based on actual logs and video pixels. `scripts/plot_research.py` generates four English PNG/PDF pairs in `figures/`; `scripts/analyze_results.py` writes aggregate and limited annotation metrics. Main logs are `logs/compare_wide.log`, `compare_target_roi.log`, `full_target_roi.log`, `quiet_confirmation.log`, `fault_target_roi.log`, `cotracker_online.log`. Raw per-frame decode/stage/received-state histories are in each `results/` run. `INPUT_OUTPUT_SHA256.json`, Git history and `evidence/execution_identity.json` identify inputs, source and artifacts. Visual checks are recorded in `evidence/visual_inspection.json`.

No video was uploaded, no camera firmware changed, no original dependency upgraded, and no instrument output was changed. The broad timing comparisons and all failed results remain locally available.

Bounded replay/GPU test wall times total **145.39 s**, excluding installation/downloads, input inspection, model loading and report/figure generation; full per-run accounting is in `evidence/test_wall_times.json`.
