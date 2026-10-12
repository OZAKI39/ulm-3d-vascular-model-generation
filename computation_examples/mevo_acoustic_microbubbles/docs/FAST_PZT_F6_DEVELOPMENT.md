# F6 — prescribed-target navigation and bounded physical-model checks

## Identity and preservation

This is an additive development from fixed F5 sync commit
`c8ecf9deece1bb0e86b0d0f1c574b984a30ebbd0`, on branch
`dev/fast-pzt-f6-target-navigation-20261010`. Remote sync branches were inspected
before work; no newer MEVO sync snapshot superseded the requested baseline.
The isolated WSL worktree is
`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f6_20261010`.
The original local project and other worktrees were left in place. All160 inherited
source files checked against Git match exactly. The13 current F5 input hashes
match the supplied manifest; old results and default F5 execution are retained.

The server was freshly inspected on2026-10-10, not inferred from migration notes:
RTX4090, driver550.120, Ryzen7950X3D, Python3.12.3, Warp1.18.0+cu12,
NumPy2.2.6, SciPy1.15.3 and h5py3.16.0 in the existing isolated compute environment.
No driver or compute-package upgrade, CFD solve, k-Wave solve, or table rebuild
was performed. The old server launcher was not used. New supervisor jobs run
from immutable per-commit F6 source directories, with one GPU science job at a
time; the two radial references used two CPU workers.

| Evidence | Actual execution source |
|---|---|
| Final admission: gradient, regression, RPY, adaptive, batch | `2ec976952edf9bad1306e90c5a82555c7c943ef9` |
| Formal17-case navigation campaign | `2caf8d8884c8b1cb2bde5dc7fb9ec5b5d0be30c2` |
| Independent low-pressure APECSS reference | `99fcbfc0bd7cadb9940e91a319af60b0a692150e` |
| Stored-table reference comparison | `788085892` (full identity in result JSON) |

Each result's `IDENTITY.json` records full source SHA-256 and the input hashes.
Postprocessing is independently versioned and never replaces the execution SHA.
See [physical specification](FAST_PZT_F6_PHYSICS_MODEL.md),
[registered configuration](../configs/fast_pzt/f6/protocol.json),
[execution/readout guide](../scripts/fast_pzt_f6/README.md), and
[historical preservation audit](../results/fast_pzt/f6/HISTORY_PRESERVATION.json).

## Implemented functions

| File in `src/mevo_mb/fast_sim/` | Behavior and qualification |
|---|---|
| `consistent_field.py` | Same complex trilinear pressure and analytical derivative on CPU/GPU; strict finite-domain rejection; explicit one-sided faces; old nodal-gradient mode retained |
| `navigation.py` | Four bounded source amplitudes and three relative phases; static optimization and25ms centroid feedback; no imposed bubble/group velocity |
| `rpy.py` | Independent equal-radius free-space RPY GPU mobility with self/overlap branches; old lubrication bypassed |
| `f6_engine.py` | Explicit gradient/mobility/correction flags, original F2 force/contact integration, unchanged contact quality gates |
| `adaptive.py` | Optional error-controlled step doubling, accepted physical-step record, diagnostic rollback, exact control-event segmentation |
| `ai_samples.py` | Group state with individual IDs/coordinates, radii, short history, shape/connectivity, controls, actual flow and quality, grouped split checks |

`scripts/fast_pzt_f6/` contains actual checks, static preflight, pilots, formal
science runner, low-pressure independent radial/translation reference, analysis,
and hash-verified archival. Existing storage datasets are preserved with new
logical HDF5 groups and hard links. Evidence is always `EXPLORATORY`.

## Frozen scientific questions and admission

Targets were saved in `preflight_v2/TARGETS_BEFORE_SCREENING.json` **before any
navigation trajectory**. The formal `REGISTERED_PROTOCOL.json` adds the feedback
rule, admitted timesteps and final split policy; it keeps the targets fixed.
Both free-space environments use the exact archived mature32-bubble F2 initial
state and a target10µm in−Y. The original Y-CFD case uses exactly the F5 frame at
300ms and a target10µm in+Z. The no-flow and slow-affine horizons are1s, the Y-CFD
horizon0.3s. The affine test is u0=(2,0,0)µm/s with dux/dy=0.05/s; it is not CFD.
Only Y-CFD enables the original open Y wall and unchanged stored water CFD.

All variants are paired at identical individual coordinates. Tests compare sound
OFF, unchanged control, phase alignment at the prescribed target, an eligible
bounded static optimizer, and eligible feedback. The source cap stays
291.237687721856Pa per channel. Static preflight rejects the Y optimizer because
it still gives negative initial target progress; no Y feedback is admitted.
This limits the conclusion to tested controls; it is not a global proof of
uncontrollability for every possible waveform.

Acceptance was fixed before formal integration: target error≤2µm, continuous
sampled residence≥50ms, connected fraction≥0.9 throughout residence, complete
finite trajectory, all32 IDs active, zero error bits, overlap≤0.01R=17.8nm,
per-pass contact correction≤0.1R=178nm, and unchanged geometric wall budgets.
Control displacement must differ from both OFF and unchanged by>0.5µm and
>three times the key timestep centroid discrepancy. Residence is sampled every
0.5ms; it is not a continuous-time proof between samples.

The key no-flow feedback is paired at10/5µs over the full horizon. Limits are
0.5µm for centroid/target-error changes,5% for Rg and0.1 for largest fraction,
with both trajectories numerically valid. Individual packed positions need not
be bitwise identical. No threshold or target was changed to make a case pass.

Initial10µs phase-aligned pilots exceeded the17.8nm overlap budget (24.36nm);
the5µs trial also failed (25.23nm later in its trajectory). Slow-feedback10/5µs
pilots reached18.57/18.37nm. One bounded2.5µs diagnosis completed30ms without
that failure, so these three controls were preregistered at2.5µs before the
formal campaign. Other cases remain10µs. These failed pilots are retained and
are not valid training samples or successful trajectories.

## Gradient and hydrodynamic evidence

Analytic multilinear pressure, in-cell finite differences, CPU/Warp agreement,
cell-face/domain checks and exact old F5 ON/OFF1ms prefixes pass. The force change
from correcting interpolation is substantial at actual F5 positions: relative
L2 changes1.743 at initial F5 and1.801 at its300ms state, with mean-direction
cosines0.813 and0.860. The mature F3-region change is0.354 (cosine0.945).
These are **model interpolation contrasts**, not measured errors against exact
acoustics. Consistency does not establish spatial grid convergence.

RPY passes self-mobility, two-body and overlapping-branch reciprocity, CPU/GPU,
symmetry and positive-eigenvalue checks for the tested small configurations.
It supplies uf+MFac without double-counting self mobility or old lubrication.
The production contact projection still uses the inherited diagonal metric.
Wall hydrodynamics, stresslets, rotational coupling and a consistent combined
RPY/lubrication/contact solve remain outside scope.

Two new APECSSv1.8 references use305 and2038.7896Pa at2.5MHz and the same Tu-fit
manifold/sigma0 assumption as the production model. ON0–48, OFF48–72, ON72–96
cycles resolves38.4µs of local switching. The single-bubble standing-gradient and
two-bubble10R far-field examples distinguish resolved-cycle mean/R0, instantaneous
Stokes/R(t), conservative added mass, and added mass plus no-slip history.
All selected signed endpoints change by<2% on128→256 samples/cycle refinement;
the largest change is0.685%. The radial solve itself uses≤1/512 cycle and1e-10
tolerance. This is a local translation discretization check, not a new radial
model validation or full-navigation history calculation.

An additional read-only audit compares the reference mean with the **actual stored
production table**: primary displacement differs by about0.28%, pair displacement
by0.48–0.64%. The history-inclusive versus fast displacement ratio nevertheless
depends strongly on pressure, interaction and switching window (roughly0.41–3.63
in these references). One ratio cannot correct a0.3–1s moving32-bubble trajectory.
Correction remains `none` / `UNCALIBRATED`; zero drive and unsupported conditions
are never extrapolated into a fitted multiplier. The old22–163% high-pressure
result is not reused as a F5/F6 error band.

## Formal outcome and complete-data gates

All17 preregistered full trajectories completed with numerical quality PASS.
The process body took1170.855s (19.514min), or52.269 valid complete trajectories
per reserved single-GPU hour for this **mixed** set of32-bubble0.3/1s horizons
and2.5/5/10µs steps. That rate includes preparation, observed compilation, control,
simulation, metrics and HDF5, but excludes interpreter startup, archival and figures.
The no-flow10µs feedback alone took49.609s for integration and51.285s including
its setup/postprocessing; its5µs counterpart took95.190s for integration.

The principal outcomes are distinct from numerical acceptance:

| Environment / control | Final target error (µm) | First arrival (ms) | Sampled residence (ms) | Navigation criterion |
|---|---:|---:|---:|---|
| No flow / OFF | 10.000 | — | 0 | Fail |
| No flow / unchanged | 7.096 | — | 0 | Fail |
| No flow / phase aligned | 17.289 | — | 0 | Fail: moves away |
| No flow / static optimized | 7.534 | — | 0 | Fail |
| No flow / feedback | 1.281 | 783.5 | 216.5 | Pass, specified model only |
| Slow affine / OFF | 10.198 | — | 0 | Fail |
| Slow affine / unchanged | 7.327 | — | 0 | Fail |
| Slow affine / phase aligned | 17.385 | — | 0 | Fail |
| Slow affine / static optimized | 7.773 | — | 0 | Fail |
| Slow affine / feedback | 1.955 | 957.5 | 42.5 | Fail: residence<50ms |
| Y-CFD / OFF | 26.599 | — | 0 | Fail |
| Y-CFD / unchanged | 27.031 | — | 0 | Fail |
| Y-CFD / phase aligned | 27.472 | — | 0 | Fail |

The no-flow feedback makes8.719µm target-direction progress while all32 bubbles
remain connected. OFF makes essentially zero progress; unchanged makes4.377µm
with an additional4.329µm transverse drift. The slow-flow feedback retains all32
bubbles but has1.463µm residual axial displacement and misses the dwell threshold.
The horizon is not extended after observing this near-miss. Y-CFD transport is
dominated by advection; the original static control maintains cohesion without
reaching the transverse target. The OFF branch ends with a largest component7/32,
whereas the two acoustic controls end32/32. No case exits the admitted domain.

The key10/5µs pair agrees to0.0009449µm in maximum common-time centroid distance,
1.8852% in Rg and0 in connected fraction; both arrive at783.5ms and satisfy all
unchanged quality budgets. This establishes resolution of the tested group-level
result, not unique microscopic packing or physical validation.

RPY changes the scientific outcome: the same feedback **policy**, still based on
a primary diagonal-Stokes drift model, never enters the target sphere, reaches
a minimum error7.542µm and ends18.987µm away. Its final displacement is approximately
(17.085,−14.364,7.041)µm while connectivity remains32/32. Control waveforms respond
to each model's own centroid, so this is a closed-loop robustness comparison,
not identical replayed waveforms or an optimally retuned RPY controller.
The old nodal-gradient variant reaches1.878µm but first arrives at961.5ms and stays
only38.5ms, also failing the frozen dwell criterion. Thus **navigation success is
not robust to these physical/discrete-model choices**. The local history reference
cannot resolve which production mobility/translation model is physically accurate.

F5 C is the only new release branch: the exact archived300ms state with old F5
controls, old gradient, sound ON and CFD ON for another300ms. It passes numerical
checks and ends32/32 connected, whereas the archived sound-OFF/flow-ON A ends7/32.
A and B were not rerun. This supports acoustic cohesion in this numerical example,
not gas fusion or a universal release law.

See the closed-data analysis JSON and English PNG/PDF figures under
`results/fast_pzt/f6/navigation_v1/analysis/`. Numerical acceptance and supported
navigation cases are separately encoded in `FINAL_GATE.json`.

## Data, efficiency and remaining scope

The AI interface preserves source, medium, flow, bubbles, control, trajectory,
morphology and quality. The source identity, actual events and physical dt vector
are retained. Rejected adaptive trials are distinct from accepted steps; central
PBD passes never advance physical time. The affine adaptive check exercises an
actual rejection/rollback; a separate17µs noncommensurate control switch is reached
exactly. Feedback plus adaptive integration is explicitly not admitted yet.

All paired no-flow/affine variants belong to one train group; the Y-target family
is test. Validation remains empty because there are only two independent target
and initial-state families. Pilot files are separate from the formal learning pool.
Do not split frames or paired controls to fabricate independent validation.

The bounded4-bubble/5ms throughput check produced identical serial/batch positions
within1e-15m. Two serial samples took0.23454s; the batch took0.11748s, including
per-run setup, metrics and HDF5 but excluding common input preparation. The
corresponding30698/61289 valid samples per GPU hour describe that tiny workload
only. Its full process body took0.70581s. They are not32-bubble/1s navigation rates.
Formal throughput includes common preparation, graph compilation as observed,
all controls, simulation, output and failed attempts, and is reported separately.
CUDA-event elapsed time includes synchronization and host-idle gaps around run;
it must not be called exclusive GPU active-kernel time.

The17 formal cases contain2,420,000 physical steps over14.2s simulated time; the
sum of CUDA-event run spans is1155.915s. The maximum post-contact overlap is17.012nm
against17.8nm, and maximum individual contact-pass correction29.997nm against178nm.
Maximum geometric wall correction is2.981nm; minimum queried wall gap is about
−2.22e-19m (roundoff). Contact corrections remain physically significant: in the
key feedback they sum to about39.65% of the unconstrained free-path length. Passing
the budgets does not validate dense packing, near-field forces or wall drag.

No evidence here is `MODEL_CHECKED` or `PHANTOM_VALIDATED`. Water parameters,
prescribed one-way CFD, pairwise near-field closure, discrete contact, unresolved
wall hydrodynamics and the uncalibrated translation law constrain physical use.
Connectivity measures a geometric graph, not bubble fusion or gas fragmentation.

Still unimplemented/unvalidated are a calibrated long-time translation correction,
a controller using collective mobility, a consistent combined RPY/lubrication/
contact solve, feedback plus adaptive integration, and an independent third
scenario family for AI validation. No AI model has been trained. The next useful
experiment is to preregister an RPY-aware control policy on the same target and
test whether it corrects the observed transverse drift. It must remain separate
from the frozen campaign and still needs physical reference evidence. Extending
the slow-flow horizon after its near-miss would not retroactively make this trial
pass. Acoustic mesh and near-wall CFD physical convergence remain unresolved.

Closed-data checks pass for all17 files, including logical groups, controls,
fixed IDs, finite coordinates/steps, short-history reading and grouped partitions.
The [final gate](../results/fast_pzt/f6/navigation_v1/analysis/FINAL_GATE.json)
separates `validation_acceptance=PASS` from `SUPPORTED_ONLY_FOR_LISTED_CASES`.
The archive verifies288 server result files (286,050,081bytes) and410 distinct
recorded source blobs against their actual Git commits. Source/build dependency
trees are omitted; no scientific trajectory is omitted. The17 full trajectories,
failed pilots, reference curves/radial records, protocol, logs and figures are
stored in WSL; this does not rely on the previously deleted old-server backup.

After acceptance, `/home/lzy/projects/microbubble_aggregation.md` received a
minimal update in its existing prelude, first act, third act and epilogue: current
input scope/gradient, low-pressure reference limits, paired navigation/cohesion,
model sensitivity, timestep/performance/data and next steps. Its original headings,
formulae and older valid records were retained. Backup is
`/home/lzy/projects/microbubble_aggregation.md.bak_F6_20261010`; before/after copies,
reviewed diff and hashes are in `results/fast_pzt/f6/progress_document/`.
No universal navigation or physical calibration completion is claimed.
