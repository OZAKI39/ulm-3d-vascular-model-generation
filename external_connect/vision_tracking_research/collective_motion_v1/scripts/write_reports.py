"""Generate reports from saved experiment records; never invent missing measurements."""
import sys,json,subprocess,importlib.metadata as md,shutil
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *

def load(name):return json.loads((ROOT/'results'/name).read_text())
def write(name,text):(ROOT/name).write_text(text.strip()+'\n')
def fmt(x):return 'NOT_AVAILABLE' if x is None else f'{x:.4g}'
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])
m=load('metrics.json');test=m['truth_tests'];piv=test['PIV_multipass'];cpd=test['CPD'];ai=load('ai_comparison.json');cam=load('camera_status.json');real=load('real_observation.json');slow=load('causal_replay.json');fast=load('causal_replay_fast.json');sim=load('simulation_bridge.json');selection=load('selection.json')
receiver=json.loads((ROOT/'evidence/camera_readonly/windows_receiver.json').read_text())['receiver']
def active(r):
    p=r['producer'];dt=p[-1]['receive']-p[0]['receive']+(p[-1]['pts']-p[-2]['pts'])
    return {'span_s':dt,'message_hz':len(r['consumer'])/dt,'motion_hz':sum(x['valid']['translation_vx'] for x in r['consumer'])/dt,'consumer_rejected_count':sum(not x['fresh'] for x in r['consumer'])}
aslow,afast=active(slow),active(fast)
dump(ROOT/'results/runtime_active_summary.json',{'multipass':aslow,'singlepass_candidate':afast,'definition':'Rates over producer active replay span including final source-frame interval; startup-inclusive rates retained in raw runtime JSON. Valid means algebraic/quality/freshness checks, not real scientific accuracy.'})
overall_rows=[]
for name in ['PIV_multipass','CPD']:
    q=test[name]['overall'];overall_rows.append([name,q['n'],fmt(q['all_six_valid_fraction']),q['false_accept_any_parameter_n'],fmt(q['false_accept_any_parameter_fraction_of_accepted']),fmt(q['global_nonaffine_correct_rejection_fraction'])])
comparison=table(['Method','Endpoint tests','All-six output fraction','Accepted failures','Failure / accepted truth outputs','Non-affine correct rejection'],overall_rows)
param_table=table(['Parameter','MAE','95th abs error','Sign fraction (nonzero truth)','Nominal 95% interval coverage'],[[k,fmt(v['mae']),fmt(v['abs_error_p95']),fmt(v['sign_correct_fraction_nonzero_truth']),fmt(v['nominal95_coverage'])] for k,v in piv['overall']['per_parameter'].items()])
window_table=table(['Frames','PTS span ms','Output fraction','MAE vx / vy (px/s)','MAE divergence / vorticity (1/s)','90% step response ms'],[[n,fmt((n-1)*.0333*1000),fmt(piv['windows'][str(n)]['all_six_valid_fraction']),'/'.join(fmt(piv['windows'][str(n)]['per_parameter'][k]['mae']) for k in NAMES[:2]),'/'.join(fmt(piv['windows'][str(n)]['per_parameter'][k]['mae']) for k in NAMES[2:4]),fmt(m['step_response'][str(n)]['first_90pct_step_response_s']*1000) if m['step_response'][str(n)]['first_90pct_step_response_s'] is not None else 'RIGHT_CENSORED at 266.4 ms'] for n in P['windows']])
matched=test['matched_AI_comparison']
ai_table=table(['Method','Pairs','All-six fraction','MAE vx / vy','Accepted false fraction'],[[k,v['n'],fmt(v['all_six_valid_fraction']),'/'.join(fmt(v['per_parameter'][j]['mae']) for j in NAMES[:2]),fmt(v['false_accept_any_parameter_fraction_of_accepted'])] for k,v in matched.items()])
runtime_table=table(['Replay pipeline','Processed / input','Software p50 / p95 / max ms','Fresh motion Hz, active span','Fresh motion Hz, startup included'],[[n,f"{r['processed_count']} / {r['input_count']}",'/'.join(fmt(r['software_processing_ms'][k]) for k in ['p50','p95','max']),fmt(a['motion_hz']),fmt(r['valid_update_rate_hz']['translation_vx'])] for n,r,a in [('PIV 2-pass',slow,aslow),('PIV candidate 1-pass',fast,afast)]])

write('README.md',f'''# F10U-Vision / Collective Motion Perception

An executed, bounded observation prototype, built on F10T-Vision-2. **No persistent bubble IDs and no hardware actions.**
Development baseline `35755720a01de6ac6d6f14ac38e2a3ca84260236`; branch `dev/f10u-collective-motion-20261017`.

Start with [FINAL_RECOMMENDATION.md](FINAL_RECOMMENDATION.md) (Chinese 23-question report), [ACTIVE_CODE_INDEX.md](ACTIVE_CODE_INDEX.md), and [FINAL_GATE.json](FINAL_GATE.json).
Five PNG/PDF figures are in `figures/`; two fully decoded MP4 files in `videos/`. All estimates and failures in `results/`. The verbose real-observation JSON is tracked losslessly as `results/real_observation.json.gz`; a fresh checkout can restore it with `gzip -dk results/real_observation.json.gz`.

{comparison}

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
''')
write('ACTIVE_CODE_INDEX.md','''# Current entry index

| Entry | Purpose |
| --- | --- |
| src/distribution.py | Read-only adapter to original OpenCV detector; center, density, covariance; no tracker |
| src/methods.py | Official OpenPIV and pycpd calls, quality masks |
| src/affine.py | Huber affine displacement fit, homogeneous matrix logarithm, rank/coverage/uncertainty |
| src/state.py | M4 messages, independent valid masks, stale/discontinuity gate and local ZMQ mock |
| src/synthetic.py | Real image patches plus explicit analytic visual motion |
| configs/protocol.json | Preregistered seeds, windows, limits and budgets |
| configs/collective_state.schema.json | Machine-readable state and future action-event schema |
| scripts/benchmark.py | Development selection and frozen independent synthetic test |
| scripts/real_observation.py | Offline causal analysis, finite camera check, latest-slot replay |
| scripts/ai_comparison.py | One official SEA-RAFT-S checkpoint, 24 pairs, no training |
| scripts/simulation_bridge.py | Read-only actual F10S HDF5 projection |
| scripts/window_response.py | Causal step-response diagnosis |
| scripts/supplement_checks.py | Local defocus, false candidates, ZMQ fault injection and schema |
| scripts/analyze.py | Truth metrics with explicit accepted/rejected denominators |
| scripts/visualize.py | Five figure groups and two PTS-preserving videos |
| scripts/write_reports.py | Reports assembled from recorded data |
| scripts/finalize.py | Integrity checks, progress-document backup/diff, output hashes |

No temporary historical entry was proven both unused and disposable, so no old code or data was deleted. Old GUI/annotation/physics entry points remain in their original index; this index is the F10U start point.
''')
write('SOURCE_AND_REUSE_AUDIT.md',f'''# Source and reuse audit

Vision repository: `{ROOT.parent}`. Baseline `35755720a01de6ac6d6f14ac38e2a3ca84260236`. Tracked baseline clean; numerous old untracked GUI evidence files preserved (`evidence/git_status_before.txt`). No applicable AGENTS.md in the workspace ancestors. No subagents used.

Read prior supervised report, frozen human masks/manifests, held-out 49/50 report, OpenCV and Norfair interfaces, state freshness/mock code, Windows bridge and original LatestFrameStream. OpenCV held-out center F1 0.862 applies only to previously reviewed ROIs, not full-frame boundaries or arbitrary experiments. No labels or detector thresholds were tuned here. No Cellpose/SAM/SAM2 training.

`src/distribution.py` calls the original `scripts/replay_benchmark.py::Detector('cc')`, threshold 10. It scopes/restores the in-memory ROI configuration. Norfair is a dependency of that legacy module, **no Tracker is instantiated**. Candidate quality is a contrast heuristic, not probability. Counts are visible candidates, not verified bubble identities.

Windows project `/mnt/f/CHRONOS_control`, WSL original `/home/lzy/projects/external_connect/CHRONOS_control`, and Tektronix project were located/read as applicable. Original Windows receiving code and `smart_online_tracking_v2/src/live_bridge.py` reused unchanged. No Tektronix module is imported by the new observation pipeline. Existing receiver made {receiver['connection_attempts']} built-in attempts within a single bounded invocation, received 0 frames. No network configuration changed.

True latest applicable simulator worktree: `{SIM}`; HEAD `31a44d3f342e5d0bca55f7f9b83f4d32b8189696`, branch `dev/fast-pzt-f10s-x-pair-onehot-trap-scope-20261015`, tracked clean. Separate public sync HEAD `c7363522ef492d09f7902c0c3382572a70beb68a` was not assumed to be the working research baseline. The older unsuffixed mevo directory is not this F10S worktree.

Input video `{VIDEO}`, SHA-256 `{sha(VIDEO)}`. 1476×1120, recorded PTS cadence 0.0333 s (~30.030 Hz); no microscope px/um or exposure timestamp calibration. Protected input hashes before/after are in `evidence/`; source and output hashes in `INPUT_OUTPUT_SHA256.json`.

Safe cleanup: **zero deletions**. No active temporary code was established to be unreferenced and without scientific value. No physics runs, hardware controls, LED tests, SAM2 or world-model training. New independent Python 3.12 environment only. Failure logs and the first OpenPIV masked-array API failure remain saved. The API fault was repaired using official masked arrays before independent testing.
''')
write('METHODS_AND_LITERATURE.md','''# Methods, conventions and primary sources

Image x points right and y down. Coordinates are ROI-local pixels. The fixed geometric ROI center xc is the reference for v0, and is output separately from the visible candidate mean. No microscope scale is available: v0 in px/s, gradients in 1/s. Positive omega means clockwise rotation on a displayed image.

`v(x) = v0 + A(x-xc)`; `delta=a11+a22`; `omega=a21-a12`; `shear1=(a11-a22)/2`; `shear2=(a12+a21)/2`. Local rigid angular velocity is omega/2. Delta is 2D divergence, not a 3D contact/aggregation rate.

Fit endpoint displacements using SciPy least_squares with Huber loss, normalized spatial design, then form H=[[F,b],[0,0,1]]. Compute L=logm(H)/actual_dt, A=L[:2,:2], v0=L[:2,2]+A*xc. Reject nonpositive time, noninvertible/ill-conditioned or complex-log transforms. This is an interval-equivalent constant generator, not instantaneous motion during acceleration. An analytic mixed finite rotation/translation test verifies the conversion; simply using (F-I)/dt would have biased it.

Require >=8 vectors, rank 3, normalized design condition <=15, hull/ROI area >=0.18, RMS residual <=1.5 px. Reject low texture, peak quality and local median outliers before fitting. Four spatial quadrants are reported; bootstrap uses quadrant resampling (24 nominal replicates) because overlapping PIV vectors are correlated. Measured coverage is below 95%; these intervals are diagnostic and **not calibrated confidence bounds**. Numerical support checks still admit some physically wrong flow fields.

Candidate means/covariances describe visible support. In general dSigma/dt=A Sigma+Sigma A^T gives only three symmetric equations for four unknown A components, and requires the same material population with negligible births/deaths/visibility changes. Isotropic covariance is unchanged by arbitrary rigid rotation. Therefore statistics-only motion values are null; no full six-parameter recovery is claimed.

Official implementations and sources:

- [OpenPIV repository](https://github.com/OpenPIV/openpiv-python), [documentation](https://openpiv.readthedocs.io/en/latest/), [examples](https://github.com/OpenPIV/openpiv-python-examples). Installed **0.25.3, GPLv3**. Actual API has no Rust/backend selector; Rust comparison NOT_RUN for this version. We used `extended_search_area_piv`, `first_pass`, `multipass_img_deform`, validation and filters. Final replaced invalid vectors are explicitly excluded, not presented as measurements. The newer repository's speed claims are not our measurements.
- [pycpd](https://github.com/siavashk/pycpd), [API](https://siavashk.github.io/pycpd/), [Coherent Point Drift paper](https://arxiv.org/abs/0905.2635). Installed **2.0.0, MIT**. Official `AffineRegistration`, source previous frame Y, target current frame X, row convention Y@B+t; our column convention F=B.T. Three initial rotations diagnose ambiguity; convergence does not establish truth.
- [SEA-RAFT official code](https://github.com/princeton-vl/SEA-RAFT), [paper](https://arxiv.org/abs/2405.14793), [official small checkpoint](https://huggingface.co/MemorySlices/Tartan-C-T-TSKH-spring540x960-S). Code and weight model card BSD-3-Clause. Only this AI model was deployed.
- [WAFT official code](https://github.com/princeton-vl/WAFT) was reviewed as an alternative; not installed or run. SEA-RAFT-S was selected for bounded two-frame inference and lower setup burden on this 8 GB GPU, not because natural-image accuracy implies microscope accuracy.
- [SciPy least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html), [matrix logarithm](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.logm.html) supply the fit and continuous-time conversion.

Dependency license texts/pins are retained under `evidence/licenses/`; third-party code and large weights are not tracked in the research Git repository.
''')
write('COLLECTIVE_STATE_SCHEMA.md','''# M4 collective observation interface

Machine-readable schema: `configs/collective_state.schema.json`. Actual message: `evidence/state_example.json`. Implementation: `src/state.py`.

The payload includes source/session/frame ID, optional capture timestamp, host receive/processing timestamps on one monotonic clock, source PTS, ROI, calibration and image quality. Real recording PTS are not camera exposure timestamps. Live receiver time, if available, is explicitly host decode time, never silently an exposure time. Camera-to-controller end-to-end latency remains NOT_MEASURED.

Distribution: visible count, center, covariance, eigenvalues, principal axis and orientation (null if near isotropic), Rg, axis ratio, density and boundary fraction. Area validity currently means geometric ROI coverage, not a validated in-focus/PDMS mask. Visible point statistics are not segmented physical bubble-group boundaries. Density bandwidth 12 px and candidate/px² units are explicit. Density/count changes can reflect ingress/egress or missed detections.

Motion: six values plus six booleans and six reasons. Invalid = null, never old value or zero. Support size, hull coverage, normalized condition number, residual, uncertainty method and interval are retained when available. CPD intervals are multi-start empirical ranges, not confidence intervals. Bootstrap absent in timed live/replay mode is explicitly null, not a fabricated precision estimate. `valid=true` means algorithmic support/quality/freshness only; `scientific_accuracy_verified=false` and `training_ready=false` remain.

Future actions: channel, actual state, amplitude and unit, phase degrees, frequency Hz, command/execution time, synchronization quality, pressure-calibration flag. Every run used `action_event=null`; there is no actuator connection. Acoustic Pa-equivalent is never assigned to an AFG voltage field.

Local ZeroMQ queues have capacity one, no unbounded backlog. Duplicate/old/PTS-reversed/stale messages are rejected. When no new image arrives, `ConsumerGate.current()` expires the previous state after 100 ms. M4 freshness invalidates the distribution as well as motion. Consumer poll watchdog records the expiration. Resets require causal warmup and do not reuse pre-discontinuity frames. `tensor_input()` returns float32 values (NaN for missing) and a separate boolean mask; it never silently fills missing values with zero. PyTorch callers can use `torch.as_tensor(values)` and `torch.as_tensor(mask)`; no model is trained.

Unit and actual local-message tests cover partial field validity, null preservation, duplicate/time reversal, stale observations, dropped frame warmup and no-new-frame expiry. `results/mock_fault_tests.json` records the transport test.
''')
write('PIV_VALIDATION.md',f'''# OpenPIV validation

Installed OpenPIV 0.25.3 (GPLv3), NumPy/SciPy backend. Rust API absent; no Rust timing claim. Single pass: 48 px windows, 24 overlap. Official two-pass: 64→32 px, 32→16 overlap; peak-to-peak threshold 1.2, original texture SD >=2, local median outlier rejection. Invalid vectors may be replaced inside the official deformation predictor but are excluded from published final vectors and fitting. Coarse flags are used internally and final invalid-vector flags are retained; final multipass numerical SNR is not exposed by the official wrapper and is not invented. Current adapter also exposes coarse-pass SNR.

Development seed 3101 compared raw, local Gaussian-background subtraction, causal median fixed-background suppression (up to 5 prior frames + current), candidate enhancement and official multipass. Selected `{selection['mode']}` using the preregistered error/availability rule. The first engineering masked-array failure was fixed and rerun on development only, preserved as `benchmark_development_api_failure.json`. Test thresholds unchanged.

{comparison}

{param_table}

These overall errors include difficult cases and only accepted truth-bearing outputs. 4 accepted local-opposite cases had no single affine truth and are counted separately; false accepted estimates were not discarded. A small residual does not prove material motion: stationary background may produce a coherent but wrong fit. Fully raw images failed the chosen correlation-quality gate; preprocessing is necessary for this recording and these images. Candidate and fixed-background versions sometimes suppress slow-moving objects too.

Clean-case mean absolute errors, six parameters in order:

{table(['Case','MAE vx, vy, delta, omega, shear1, shear2'],[[c,', '.join(fmt(v['mae']) for v in piv['cases'][c]['per_parameter'].values())] for c in ['translation','rotation_cw','rotation_ccw','contraction','expansion','stretch','shear']])}

Fixed PDMS-like retained background, dropout, sparse and locally opposing motions reveal failures. Sparse cases are not universally mathematically unobservable: enough accidental window texture can still pass support tests, which is a limitation, not a reason to mark them successful. The supplementary local-defocus/distractor checks are post-benchmark diagnostics and did not change gates. Nominal interval coverage fails calibration; do not feed these values into a controller as certified measurements.
''')
write('CPD_VALIDATION.md',f'''# CPD independent comparison

Official pycpd 2.0.0 AffineRegistration; 60 iterations, tolerance 1e-4, w=0.1; center/scale-normalized points; initial rotations 0, -0.18 and +0.18 rad. Previous point set is source Y; current is target X. Transform converted consistently to the same six parameters. Independent implementation source convention is checked against a known affine map.

Require >=8 detected points, count change <=25%, convergence, residual <=3 px, axis ratio >=1.15 and nondegenerate design. Multi-start parameter spread above preregistered tolerances invalidates the corresponding field. Isotropy is a conservative warning, not a proof that every isotropic-looking irregular point set is ambiguous. Empirical ranges are not calibrated confidence intervals.

{comparison}

CPD sometimes has lower errors on the subset it accepts, and rejects more cases. In the matched AI subset its 10 accepted cases have no tolerance failures; that does not prove it dominates PIV on all 18 inputs. On the full test it has 12 erroneous accepted truth estimates and admits 6/12 non-affine opposing-motion cases. Static candidates can produce a wrong correspondence even when point counts stay constant; high dropout often breaks convergence/residual tests. Small nearly symmetric sets invalidate the material rotation/strain interpretation. This is an independent diagnostic, not a mandatory real-time ID tracker.

Full convergence, residual, count, start sensitivity, time and null reasons: `results/benchmark_test.json`. Additional local defocus and extra candidates: `results/supplement_stress.json`. No persistent IDs are used.
''')
write('AI_FLOW_COMPARISON.md',f'''# Bounded AI flow comparison

Selected official SEA-RAFT-S, commit `{ai['official_git']}`, model `{ai['model_id']}`. Official weights SHA-256: `{ai['model_files_sha256']}`. Code/model BSD-3-Clause. WAFT reviewed only; no second AI flow implementation deployed. Natural-image pretrained weights were not trained on these microscope images.

Actual environment {ai['torch']}, CUDA {ai['cuda']}, {ai['device']}. Independent copied CUDA dependency files, no modification to old annotation/Windows environments. Model consumes real two-frame inputs. 256×256 synthetic and 256×384 true recorded ROIs, gray repeated across RGB channels, native resolution, 4 iterations; official normalization/padding. One warmup, forward + reverse flow (forward/backward error <1.5 px), original texture rejection, same Huber/log-transform fit as PIV.

24 pairs: 18 exact synthetic endpoints (six cases × 3 windows, seed 7201) + six recorded pairs. All methods below use the exact same 18 synthetic input pairs and truth. Different acceptance subsets remain explicit.

{ai_table}

Actual final profiled p50/p95/max forward+reverse pipeline: {ai['timing_ms']['p50']:.2f}/{ai['timing_ms']['p95']:.2f}/{ai['timing_ms']['max']:.2f} ms. Fit is recorded separately, not hidden inside these numbers. Peak PyTorch allocated/reserved {ai['allocated_peak_mb']:.2f}/{ai['reserved_peak_mb']:.2f} MiB; process RSS and device-wide GPU utilization in `resources` (desktop usage is not isolated). GPU context/desktop memory is larger than the tensor allocation. Initial run preserved; exact-input repeat added process/device resource measurements, not retuning. Both finite runs are retained; no large pressure test.

**No accuracy advantage established.** AI accepts 17/18 synthetic pairs but 6/17 accepted estimates exceed a preregistered error tolerance. PIV accepts 14/18 with 0/14 such failures in this subset. AI can follow fixed background or average tiny structures. Its smoother field is not evidence of accuracy. Retain traditional PIV as the main line; no fine-tuning or world-model training this round. Real-image AI disagreement is diagnostic only because absolute six-parameter truth is unavailable.
''')
write('WINDOW_3_5_10_BENCHMARK.md',f'''# Causal 3/5/10-frame windows

Endpoints span N-1 actual PTS intervals. Recorded cadence is ~30.030 Hz, not assumed 60 Hz. Synthetic time uses explicitly registered 0.0333 s intervals. All frames available to the method are at or before the current endpoint. No future frames or bidirectional temporal smoothing.

{window_table}

Errors condition on accepted outputs; lower error when many cases are rejected is not by itself an improvement. Detailed bias, SD, signs, CI coverage, false accepts, availability and processing times for each method/window are in `results/metrics.json`. Figure 4 shows the accuracy/availability tradeoff and a causal sudden translation step.

The 10-frame window has not reached 90% by the short test horizon, 266.4 ms after onset; this is right-censoring, not a measured infinite lag. The interval-equivalent velocity mixes old and new motion as it should. Window midpoint ages are 33.3/66.6/149.85 ms, in addition to software processing delay. Longer windows can also lose correspondence or average different local motions. No adaptive-window policy was tuned on the final test.
''')
write('SYNTHETIC_GROUND_TRUTH.md','''# Synthetic known-motion ground truth

This is visual/analytic synthesis, not an acoustic experiment. Real Chronos dark bubble image patches are placed on a flat or retained real image background. Each center moves by an exact homogeneous exponential; individual small image patches retain their appearance, so this is advected visual texture, not a deforming acoustic bubble model. Patches are from OpenCV visible candidates, not asserted physically verified bubbles.

Development: seed 3101, source frame 90. Independent synthetic test: seeds 7201/7202, source frame 840. Background/patch source differs between development and test, but only one recorded session; no cross-session generalization claim. No human test masks or consumed 49/50 labels are used for tuning. 18 cases, 16 frames each; final endpoints 9/15, windows 3/5/10. Thus 216 endpoint comparisons per method, **not 216 independent experiments**; adjacent windows and a shared background are correlated.

Cases: XY translation, isotropic expansion/contraction, clockwise/counterclockwise rotation, stretch, shear, mixed motion, fixed background, density 100 versus nominal 40, four sparse targets, severe blur, 45% random missed visibility, boundary crossing, symmetric ring, opposing local translations, brightness/noise and sudden translation. The local-opposite case has no unique single global affine truth. Static real texture may contain other apparent dark candidates and is deliberately retained as a confounder. Known underlying centers and visible generator IDs are saved solely for validation, never supplied to the observation algorithm. The detector is run independently on both endpoint images.

Saved per-case NPZ: actual images, H matrices, true centers; JSON: seed, source PTS, dt, generated visibility, case, global affine truth flag; manifests include SHA-256. Raw motion truth is computed with the same definition but directly from the saved analytic H, never fitted to the observations. Scoring uses fixed absolute tolerances [5,5,0.2,0.2,0.1,0.1] in [px/s,px/s,1/s,1/s,1/s,1/s]. They are development criteria, not experimental calibration.

`supplement_stress.json` adds explicitly local half-ROI defocus and extra false candidates after the primary benchmark, with no retuning. Video-only seed 9001, 48 frames × 6 cases gives a readable 9.5904 s sequence at original synthetic cadence; these additional calculations are labeled visualization and excluded from the held-out summary. Long synthetic motion is not presented as a physical simulator trajectory. No successful/failed case was removed.
''')
write('REAL_CHRONOS_OBSERVATION.md',f'''# Actual Chronos recording and live availability

Saved video SHA-256 `{sha(VIDEO)}`. Frame IDs 360–479, original PTS {real['rows'][0]['pts']:.4f}–{real['rows'][-1]['pts']:.4f} s. ROI [650,300,256,384]. 120 frames decoded and analyzed causally with the original detector; no persistent bubble IDs. The mask/shape/ID accuracy limitations of F10T still apply. Visible density and covariance do not measure 3D concentration or a fused physical cluster.

Real absolute six-parameter ground truth is **NOT_AVAILABLE**. The displayed original images and arrows are qualitatively consistent with a dominant image-downward motion; this visual check does not certify local divergence, rotation or strain. All 3/5/10-window numerical fits on this selected ROI pass their support criteria, which does not validate their absolute accuracy. Changes of the visible mean are irregular due to changing count, edge visibility and detection. CPD and AI comparisons are recorded, not treated as independent truth. Real out-of-focus ground truth is not annotated; focus score is diagnostic.

The existing Windows receiver was invoked once for a finite ~8 s read-only check. It made {receiver['connection_attempts']} internal attempts, received {receiver['received_frames']} frames, and reported `{receiver['last_error']}`. Network configuration and original receiver source unchanged. `LIVE_CAMERA_PASS=false`. All subsequent runtime claims are `CAUSAL_REPLAY_TEST`, not live-camera performance.

The MP4 preserves original frame PTS relative to the first displayed frame, including all 120 original frames. It shows causal endpoint estimates computed offline and explicitly says offline rendering. It is not a recording of the slower wall-clock receiver pipeline, and does not hide that runtime skipped frames. Runtime rows separately record overwritten/skipped input, stale/null observations and software latency.
''')
write('SIMULATION_OBSERVATION_BRIDGE.md',f'''# Actual saved simulator observation bridge

File `{sim['file']}`; SHA-256 `{sim['sha256']}`; F10S source HEAD `{sim['source_git']}`. Read-only HDF5 contains 3001 time/position samples for 32 fixed particle IDs across **0.003 s**, plus force/pressure samples and solver diagnostics. No explicit velocity dataset: adapter derives interval finite-difference velocity from stored positions and times, and labels it accordingly.

Original numerical basic/group position/group shape/topology/ID gates PASS for this saved run. Physical model remains `PHYSICS_UNCALIBRATED`, `HOMOGENEOUS_WATER_CONTROL_PROXY`, `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED`, `training_ready=false`. Numerical pass is not experimental validation. No new physics/GPU trajectory executed.

XY projection gives COM displacement {sim['com_displacement_um'][0]:.9f} um in X over 3 ms; fitted reference-center translation {sim['six_parameters']['translation_vx']:.6f} um/s. Full projected covariance, affine residual and upper/lower regional fits are in `results/simulation_bridge.json`. IDs are used only to establish saved simulation truth; vision does not receive them. Lower/upper residual diagnostics do not identify a need for arbitrary complex shape control from this tiny, nearly uniform signal.

Projection loses out-of-plane motion. There is no measured PSF, pixel scale, illumination/depth/visibility model for this device, so simulator-to-microscope rendering and learned adaptation are **NOT_RUN**. Do not stretch/repeat this 3 ms trajectory into seconds, and do not call old Pa-equivalent proxy source amplitudes AFG volts. Future upgrades: measured image scale/PSF/depth, optical forward model with occlusion and fixed artifacts, independently labeled camera data, action/clock schema calibration. No need for new expensive PDMS/CFD/contact physics to complete this observation adapter.
''')
write('PERFORMANCE_AND_ADMISSION.md',f'''# Measured runtime and admission

{runtime_table}

Active acquisition spans {aslow['span_s']:.4f}/{afast['span_s']:.4f} s. Startup-inclusive raw rates include decoder preroll; both denominators are reported explicitly. M4 message rates during active replay are {aslow['message_hz']:.3f}/{afast['message_hz']:.3f} Hz, while fresh six-parameter rates are {aslow['motion_hz']:.3f}/{afast['motion_hz']:.3f} Hz. Candidate distributions can remain available during motion warmup. Skipped frames: {slow['input_count']-slow['processed_count']}/{fast['input_count']-fast['processed_count']}. Consumer rejection includes both expired processing and frame-discontinuity resets, not only stale processing. Latest slots/queues are bounded at one; there is no accumulated FIFO lag.

The 2-pass p50 includes cheap warmup/reset iterations, **not** a 6 ms full multipass solve. Detailed PIV and fit-only timing is in `causal_replay.json`. Single-pass is a secondary latency probe selected from already tested development options, not a new independently tuned accuracy winner. It did not undergo a new independent held-out accuracy trial after this speed probe. No claim of 30/60 Hz complete real-time six-parameter processing.

Transport p95 {slow['transfer_ms']['p95']:.3f}/{fast['transfer_ms']['p95']:.3f} ms; CPU core equivalents {slow['cpu_core_equivalents']:.3f}/{fast['cpu_core_equivalents']:.3f}, process RSS samples retained. PIV/CPD run on CPU. SEA-RAFT tensor peak {ai['allocated_peak_mb']:.2f} MiB, measured device-wide utilization and process RSS retained; includes other desktop GPU users. AI times are forward+reverse flow plus mask/transfer; not a single-forward-only benchmark.

Measured independent synthetic test processing {load('benchmark_test.json')['elapsed_s']:.3f} s, real 120-frame all-window offline analysis {real['elapsed_s']:.3f} s. Benchmark timers exclude synthetic generation/setup; AI resource interval {ai['resources']['elapsed_s']:.3f} s excludes weight loading/warmup. No invented total GPU-hour/whole-development estimate. Install/plot/video time is not online performance.

Admission: software/math/interface PASS with explicit limits; clear-case synthetic recovery LIMITED_PASS; difficult-case false acceptance and calibrated uncertainty FAIL; real absolute motion accuracy NOT_VERIFIED; live camera unavailable; hardware/PZT/AFG NOT_RUN; CAMERA_TO_CONTROLLER_END_TO_END_LATENCY NOT_MEASURED. `training_ready=false`. Freshness does not establish scientific accuracy.
''')

answers=[
('是否复用现有代码？','已复用原 OpenCV Detector、WindowsLatestBridge 与原 Chronos LatestFrameStream，未改相机或 AFG 工程。'),
('是否取消长期 ID 依赖？','是。没有实例化 Norfair Tracker；候选点、PIV 和 CPD 均不需要单泡长期身份。'),
('中心与分布能否连续输出？','120 帧真实录像已输出候选中心、数量、密度、协方差、Rg、长短比和方向；近各向同性方向为空。它们是可见二维候选分布，不能冒充真实团边界。'),
('PIV 能恢复已知运动吗？','清晰平移、顺逆旋转、收缩、扩张、伸展、剪切测试可恢复；下表保留复杂背景/漏检的失败，不能宣称所有六参数已普遍验证。'),
('CPD 表现如何？',f"全测试全六参数输出率 {cpd['overall']['all_six_valid_fraction']:.1%}；有真值的已接受输出中 {cpd['overall']['false_accept_any_parameter_fraction_of_accepted']:.1%} 超过至少一项预设容差。清晰且非对称条件下精度较好，拒判更多，先作低频对照。"),
('AI 光流有收益吗？','实际运行官方 SEA-RAFT-S。相同 18 对合成输入上，AI 接受 17 对，其中 6 对超容差；PIV 接受 14 对，其中 0 对超容差。暂未显示精度收益，不替换传统主线。'),
('哪些参数容易失败？','固定纹理/漏检下平移和涡量会明显偏差；散度、涡量的经验区间覆盖尤其不足。不能把低残差或平滑曲线当可信。'),
('对称、漏检、进出视野影响？','对称点集使材料映射/旋转有歧义；漏检和静态假点会改变配准；边界进出使可见均值移动但未必是材料平移。三类均有实际诊断和空值原因。'),
('3/5/10 帧如何？','实际跨度 66.6/133.2/299.7 ms，PIV 全六参数输出率 91.7%/62.5%/55.6%；条件误差见下表。阶跃90%响应为66.6/133.2 ms，10帧在266.4 ms观察期内未到90%。长窗并非总是更准。'),
('哪种适合未来实时控制？','研究主线保留 PIV；候选增强单遍优先做低延迟原型，两遍作离线/低频复核，CPD作独立诊断。目前所有方法均未获得真实控制准入。'),
('哪些通过真值测试？','矩阵对数转换、坐标符号、CPD矩阵约定、清晰解析仿射运动、空值/过期/时间断裂接口已验证；困难场景和95%置信区间没有通过通用精度验收。'),
('哪些只有定性真实证据？','真实录像的运动方向和局部六参数趋势只有图像/算法间对照，缺少绝对六参数人工真值；未做PZT ON/OFF因果试验。'),
('GPU仿真提供什么？','读取真实 F10S HDF5：32泡、3001个位置/时间样本、完整3ms和已有数值门槛。速度由位置差分；X位移0.197999447 um。未新跑仿真。'),
('哪些仿真升级留以后？','显微标尺、PSF/景深/遮挡、固定纹理、真实动作时序的观测校准；没有必要在本轮重启高成本声学或PDMS开发。'),
('直播成功了吗？','没有。原Windows接收器本次0帧，LIVE_CAMERA_PASS=false；只执行一次有界只读会话，保留内部重连失败记录。'),
('更新率和p95是多少？',f"真实录像因果回放：两遍PIV p95 {slow['software_processing_ms']['p95']:.2f} ms，新鲜六参数0 Hz；单遍候选增强 p95 {fast['software_processing_ms']['p95']:.2f} ms，采集活动时段新鲜六参数 {afast['motion_hz']:.2f} Hz（含启动 {fast['valid_update_rate_hz']['translation_vx']:.2f} Hz）。不等于直播/端到端性能。"),
('哪些NOT_OBSERVABLE？','低支持/覆盖、退化设计、时间异常、过期、相关失败、CPD点数变化/对称/多起点歧义分别拒判；点集统计不唯一识别完整六参数，全部运动字段为空。仍存在错误接受并已计数。'),
('统一Schema是否输出？','已输出M4 JSON Schema、实际消息、按字段valid/reason、NaN+mask张量适配、有限ZMQ Mock和无新帧过期检查。动作仅Schema，实际action_event=null。'),
('是否执行PZT控制？','没有。未调用AFG SCPI、功放、LED、SAM2或世界模型训练。'),
('图和视频在哪里？',str(ROOT/'figures')+'（五组PNG/PDF）；'+str(ROOT/'videos')+'（两段MP4）。'),
('视频检查了吗？','两段逐帧完整解码：120/288帧；PTS严格单调，误差约机器精度，无插帧；已查看头/中/尾与各典型案例截图。'),
('哪些未完成，为什么？','真实直播因相机无画面未完成；曝光到控制器延迟、跨实验日精度、真实六参数真值、可信区间校准及PZT响应尚未验证。Rust后端不在已安装OpenPIV版本中；WAFT和训练按范围未运行。'),
('下一轮人工PZT响应验证前缺什么？','先恢复相机并校准像素尺度/时间；保留独立显微运动真值和真实固定背景反例，修正错误接受与区间覆盖；优化ROI/调度并测真实接流延迟。之后另行开展人工低风险动作记录，不能由本轮软件PASS自动获得硬件控制权限。')]
write('FINAL_RECOMMENDATION.md',f'''# F10U-Vision 实际交付与23项回答

已实现并实际运行“无长期单泡ID”的群体观测原型。保留 PIV 主线、CPD 对照，SEA-RAFT 没有显示显微场景收益。**这是完成有限研究和接口交付，不是已通过真实自动控制验收。**

{chr(10).join(f'{i}. **{q}** {a}\n' for i,(q,a) in enumerate(answers,1))}

## 完整测试与分母

{comparison}

{window_table}

{runtime_table}

原始失败、条件误差、符号、覆盖率、未观测原因均保留于results。216个端点对比来自两个测试种子、18场景和重复时间窗，不能称216个独立实验。固定背景可能产生低残差假运动；12个局部相反运动病例里PIV仍接受4个、CPD仍接受6个，说明拒判机制仍需改进。采集活动时段与含启动耗时的更新率均报告，未隐去启动与掉帧。

实现与执行基线：`3575572` → 新分支 `dev/f10u-collective-motion-20261017`。各运行JSON保存执行SHA和源码SHA；最终提交见`evidence/delivery_git.json`。输入/输出哈希、图像检查、完整解码和进度文档before/after/diff均有独立文件。未推送GitHub，因为本轮任务只要求本地开发交付。
''')
gate={'stage':'F10U-Vision','SOFTWARE_IMPLEMENTED':True,'KNOWN_AFFINE_MATH_TEST':'PASS','SYNTHETIC_CLEAR_CASE_RECOVERY':'LIMITED_PASS','DIFFICULT_CASE_GENERALIZATION':'FAIL_FALSE_ACCEPTS_RETAINED','CALIBRATED_95_PERCENT_INTERVALS':'FAIL_UNDERCOVERAGE','REAL_SIX_PARAMETER_ACCURACY':'NOT_VERIFIED','LIVE_CAMERA_PASS':False,'CAUSAL_REPLAY_TEST':True,'MULTIPASS_FRESH_MOTION_HZ':aslow['motion_hz'],'SINGLEPASS_FRESH_MOTION_HZ_ACTIVE':afast['motion_hz'],'CAMERA_TO_CONTROLLER_END_TO_END_LATENCY':'NOT_MEASURED','AI_MODEL_RUN':'SEA_RAFT_S','AI_ACCURACY_ADVANTAGE':'NOT_DEMONSTRATED','PERSISTENT_BUBBLE_IDS_REQUIRED':False,'PZT_CONTROL_EXECUTED':False,'AFG_CONTROL_EXECUTED':False,'WORLD_MODEL_TRAINED':False,'TRAINING_READY':False,'training_ready':False,'physical_3d_aggregation_inferred':False,'scientific_control_admitted':False,'figure_groups':5,'videos':2,'limitations':['Single recorded session','Visible candidate distribution, not segmented physical cluster','No micron/time exposure calibration','Algorithmic observability does not guarantee true material motion','No calibrated confidence intervals']}
dump(ROOT/'FINAL_GATE.json',gate)
license_dir=ROOT/'evidence/licenses';license_dir.mkdir(exist_ok=True)
for name in ['openpiv','pycpd']:
    d=md.distribution(name)
    for f in d.files:
        if 'LICENSE' in str(f):shutil.copyfile(d.locate_file(f),license_dir/(name+'_LICENSE.txt'))
shutil.copyfile(ROOT/'third_party/SEA-RAFT/LICENSE',license_dir/'SEA_RAFT_LICENSE.txt')
dump(ROOT/'evidence/official_sources.json',{'SEA_RAFT_git':ai['official_git'],'SEA_RAFT_url':'https://github.com/princeton-vl/SEA-RAFT','SEA_RAFT_weights':ai['model_id'],'versions':{k:md.version(k) for k in ['openpiv','pycpd','numpy','scipy','torch','torchvision','av']}})
print('required reports written; final integrity/progress checks remain')
