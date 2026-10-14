# Current acoustic–microbubble project — F10Q, 2026-10-14

Latest local scientific source, including the completed parameter/force summary:
`41ae8ed2b88f2cbbed834fffccaa0b5de6e91589`.
Source worktree:
`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10q_20261014/computation_examples/mevo_acoustic_microbubbles`.
Curated branch: `sync/mevo-acoustic-f10q-20261014`.

This supersedes the F10K handoff for current development. F1–F10P code and
conclusions remain historical records. The scientific worktree was clean at
export; this task changes packaging only.

## Model and entrypoints

Four independently controlled 7 mm PZT discs at X± and Y± faces of a
32.5 × 9 × 9 mm apparatus, a central 1.6 mm channel, two 1.5 mm top openings
and closed X ends. Motion calculations use water-filled accessible cavities,
strictly zero background flow, homogeneous-water disc acoustics at 2.5 MHz,
APECSS response and 32 individually tracked bubbles of 1.78 µm radius. The
connected initial state is prescribed, not formed by this acoustic experiment.
Source amplitude means rho*c*normal-velocity Pa equivalent, not measured
experimental PZT face pressure or voltage.

- [Active code map](F10Q_ACTIVE_CODE_MAP.md) and
  [Chinese F10Q report](FAST_PZT_F10Q_DEVELOPMENT.md).
- [True-gap audit](FAST_PZT_F10Q_TRUE_GAP_AUDIT.md),
  [sequential QP](FAST_PZT_F10Q_SEQUENTIAL_QP.md),
  [group comparison](FAST_PZT_F10Q_GROUP_COMPARISON.md) and
  [XY navigation](FAST_PZT_F10Q_XY_NAVIGATION.md).
- Protocol: `configs/straight_tube_exp1/f10q/`; opt-in solver:
  `src/mevo_mb/fast_sim/sequential_hard_sphere.py` and `water_sequential.py`.
- Protected original inputs: `results/straight_tube_exp1/f10n/`, including
  initial states, APECSS response and both original acoustic caches.

## Actual outcomes

X control reaches the pre-registered 0.2 µm target over 2 ms in this simplified
simulation: displacement 0.200896 µm and endpoint error 0.000896 µm. Position,
shape, topology and fixed-ID numerical gates pass. It is nearly translation,
not demonstrated arbitrary shape control.

Y displacement is 0.051776 µm toward its independent 0.2 µm target; endpoint
error remains 0.148224 µm. Position/topology pass, but shape and fixed-ID
convergence fail. OFF remains stationary. The old strong control advances
validly to 0.883 ms and rejects the 0.884 ms proposal with an active-set cycle.
Its full 2 ms case is not admitted; continuation and 5/10/30 ms runs were not run.
Recorded numerical time is 139.215 s, or 141.020 s including read-only numerical
verification. Hybrid solver time is not pure GPU time.

All results remain **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY /
PHYSICS_UNCALIBRATED / training_ready=false**. Contact-model dependence and
experimental calibration of the new controls are unverified. Geometric graph
connectivity does not mean physical fusion.

## Latest evidence

- `results/straight_tube_exp1/f10q/analysis/FINAL_GATE.json`: layered gates.
- `results/straight_tube_exp1/f10q/`: original trajectories, failures, inputs,
  logs, source execution archives, manifests and four final PNG/PDF groups.
- `results/straight_tube_exp1/f10q_parameter_summary_20261014/`: PZT layout,
  source amplitudes/phases, actual pressure, per-bubble mean primary force and
  maximum pair attraction over all accepted force-evaluation steps. No new
  trajectory was run. Y statistics retain the shape-convergence limitation.
- [Chinese research progress](MICROBUBBLE_AGGREGATION_PROGRESS_20261014.md):
  unchanged live-document copy; original before/after/diff audit also retained.

## Environment and reproduction boundary

Existing remote activation:
`source /workspace/mevo_acoustic_microbubbles/activate.sh`.
Dependency declarations and recorded environment/execution evidence are retained;
environments and private keys are excluded. F10Q server data was returned and
hash-verified before this synchronization. No server jobs are started or stopped.

Run the packaging verifier from the project directory:
`python scripts/verify_sync_20261014.py`. Scientific reruns need the numerical
environment and a fresh output location; do not overwrite frozen results or
treat historical RUNNING records as instructions to resume jobs.
See [scope, omissions and recovery](../sync/20261014/README.md).
