# MEVO acoustic, microbubble and vision synchronization — 2026-10-17

This is a curated **new branch**, `sync/mevo-acoustic-vision-f10u-20261017`, based on
published F10S sync `c7363522ef492d09f7902c0c3382572a70beb68a`.
It adds the completed F10T/F10U vision work and associated device interfaces.
No solver, camera, PZT, AFG, annotation GUI or training service was started for
this synchronization. Existing scientific outputs and development worktrees are
preserved; the main branch and other repository projects are not changed.

## Start here

| Component | Repository path | Source identity |
| --- | --- | --- |
| Acoustic field, four-PZT control, full-particle microbubble solver through F10S | `computation_examples/mevo_acoustic_microbubbles/` | `31a44d3f342e5d0bca55f7f9b83f4d32b8189696` |
| Chronos research, annotation, detection and F10T held-out evaluation | `external_connect/vision_tracking_research/smart_online_tracking_v2/` | Included in vision delivery below |
| No-ID collective motion, six image-plane parameters and observation interface | `external_connect/vision_tracking_research/collective_motion_v1/` | `0561769933dfa799dee0778f638a8c26270c3b95` |
| Original camera receiver and acquisition interfaces | `external_connect/CHRONOS_control/` | File SHA-256 in `COPIED_FILES.json` |
| Existing AFG interfaces and hardware test evidence | `external_connect/Tektronix_control/` | File SHA-256 in `COPIED_FILES.json`; not operated by this sync |
| Current Chinese research progress | `computation_examples/mevo_acoustic_microbubbles/docs/MICROBUBBLE_AGGREGATION_PROGRESS_20261017.md` | Exact copy; original progress document unchanged |

Use `docs/F10S_ACTIVE_CODE_MAP.md` in the simulation project and
`ACTIVE_CODE_INDEX.md`, `FINAL_RECOMMENDATION.md`, `FINAL_GATE.json` in the F10U
directory. The annotation instructions remain at
`smart_online_tracking_v2/ANNOTATION_READY.md`. Latest human-label evaluation is
in `smart_online_tracking_v2/results/heldout_49_50_20261016_v1/REPORT.md`.

## Included and excluded

All 1,098 vision Git-tracked files are retained byte-for-byte, along with
previously ignored scientific evidence: saved annotation revisions, the 11
annotated candidate images, all three frozen annotation snapshots (including
image/mask TIFF exports), result tables, success/failure logs, synthetic ground
truth, actual flow arrays, selected real frames, final PNG/PDF figures and six
small derived MP4 deliverables. The first five partially annotated frames stay
`DRAFT_UNREVIEWED`; unmarked pixels are not relabeled as audited background.

The new copied payload is approximately 220 MB before Git packing; every new
file is below 10 MB. Forty-one verbose text results are losslessly gzip-compressed.
The additional manifests describe each file's origin, size and SHA-256. Existing
F10S physics files are inherited unchanged, including convergence trajectories,
field caches, response tables and original failures. A read-only server check
confirmed all 75 archived F10S result files match their remote originals.

Excluded: raw source movies, unannotated candidate images, repeated frame dumps,
temporary segmentation proposals/preview images, model weights, Python
environments, downloaded dependency trees, VLC/TekVISA binaries, caches and Git
databases. `EXCLUDED_FILES.json` identifies individually reviewed omissions and
their hashes. `EXCLUDED_DIRECTORIES.json` lists pruned dependency/environment
directories; their sizes and every member's hash were **not** inventoried.
Omitted local originals have not been deleted. Historical source Git commits are
provenance identities, not a promise that those private development commits are
reachable on GitHub.

## Verification and restoration

From the repository root, with NumPy installed, verify the package before
restoring omitted/expanded files:

```bash
python computation_examples/mevo_acoustic_microbubbles/sync/20261017/verify.py
```

`PACKAGE_MANIFEST.json` covers the newly published files except itself and
`VERIFICATION.json` (to avoid a self-hash cycle). The parent commit identifies
the inherited tree. `COPIED_FILES.json` links exported bytes to original source
bytes; compressed files have both hashes. `VERIFICATION.json`, test logs and
`REMOTE_READONLY_AUDIT.json` record the actual checks.

To expand the 41 compressed text results without overwriting modified files:

```bash
python computation_examples/mevo_acoustic_microbubbles/sync/20261017/verify.py --restore-compressed
gzip -dk external_connect/vision_tracking_research/collective_motion_v1/results/real_observation.json.gz
```

The F10U JSON above was already published locally in compressed form, separately
from this sync's 41 conversions. Expanded bytes match their original SHA-256.

Create independent environments from each module's supplied requirements/lock
files; environments and model downloads are not bundled. F10U's frozen torch
build uses the CUDA 12.8 wheel index. Read its `AI_FLOW_COMPARISON.md` and
`results/ai_comparison.json` for the official SEA-RAFT code/checkpoint identities;
the AI comparison needs that separately restored dependency and checkpoint.
Model weights used by the annotation workbench are also separate downloads or
local training artifacts, not silently replaced by another model.

Original source code and records deliberately retain their original absolute
paths. On another machine, configure paths in a separate execution copy:
`collective_motion_v1/src/common.py` defines the source recording and simulator
paths; camera bridge code refers to the user's `F:\CHRONOS_control` receiver;
historical configurations and frozen records retain their original paths.
The shared `external_connect/` layout is preserved so relative imports resolve.
Restore the source movie with the recorded SHA-256 for full replay, and the
remaining candidate frames for complete annotation browsing. The public package
already includes the reviewed frames and frozen masks needed to inspect the
evaluation. It is not a turnkey hardware installation or an environment backup.
Run new experiments in separate output directories; do not rerun historical
delivery/finalization scripts over this evidence snapshot.

## Scientific status preserved

F10S remains `EXPLORATORY`, `HOMOGENEOUS_WATER_CONTROL_PROXY`,
`PHYSICS_UNCALIBRATED`, `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED`,
`training_ready=false`. Its limited X navigation success does not establish
arbitrary shape control or a calibrated PDMS device model.

F10U implements an observation prototype. Difficult visual cases retain false
accepts and undercovered confidence intervals. Real six-parameter accuracy is
`NOT_VERIFIED`, `LIVE_CAMERA_PASS=false`, `scientific_control_admitted=false`,
`training_ready=false`. Historical F10T live tests are separate records; they do
not turn the later F10U disconnected-camera test into a pass.
