# Fast PZT-controlled aggregation F1 — development and measured results

The standalone Warp prototype runs on the RTX 4090, forms spatial clusters from
dispersed bubbles and responds to a change in real precomputed PZT basis controls.
The delivered physics scope is **EXPLORATORY**, not quantitatively calibrated.
The machine-readable capability gates are in
[`FINAL_GATE.json`](../results/fast_pzt/f1/FINAL_GATE.json).

Base: reviewed R3 archive `5c0a0f10b19b376d9321510149cb6db50bcf4331`.
Branch: `dev/fast-pzt-aggregation-f1-20261007`.
Final scientific runs executed immutable source `9829156893339efff713483f3125848b84adbc98`;
the matched OFF/two-bubble controls used `3721bebd597164262691af766f723d245e043116`
(identical simulator source, additional control-scene runner). Every run stores source
hashes and input hashes. The original Stage 4/5 code, data, failed records and gates
were preserved. The earlier R3 worktree's untracked `SOURCE_COMMIT_R3.txt` was untouched.

## 1. Independent GPU implementation

`src/mevo_mb/fast_sim` is explicitly selected by the new scripts; no default backend
was replaced. Warp 1.18.0 executes float64 kernels on the RTX 4090, with persistent
positions, eight volume harmonics, force/velocity buffers, adjacency state and quality
counters. Full-pair forces use every unordered pair; 200 bubbles have 19,900 pairs.
Bounded trilinear field query, response evaluation and Primary force are fused.
CUDA graphs launch chunks of 32 actual steps. Control events split steps exactly.
No per-particle Python stepping, global QP, wall geometry or CFD is involved.

The isolated Python 3.12.3 environment and exact dependency lock are archived.
Driver 595.84 was retained; Warp reports toolkit 13.4 / driver CUDA 13.2 and compiled
successfully. Final-model cold module compilation took 2.295 s at the default block
size and 2.266 s for the second tested layout. Cached compilation/load time and graph
capture are recorded separately. Twelve small GPU/NumPy and invariant checks passed:
field/phase/response/Primary, harmonic pair force, OFF, force signs/reciprocity,
normal mobility, separating contact, noncommensurate control event, domain rejection,
near-field continuity, thread layout, independent batch members and component balance.
The signed repulsion check uses synthetic counterphase harmonic responses; identical
bubbles in the uniform production scene attract. It is a formula test, not a claim
that every production pressure produces repulsion.

## 2. Actual aggregation and its control

Scene A: the two-bubble ON case closes an 8 µm separation to approximately the
3.56 µm diameter; the matched OFF pair remains fixed. Scene B starts with 32 random,
nonoverlapping bubbles (seed 20261007, cloud radius 20 µm, zero flow, 1 kPa).
Secondary ON changes 32 components into one 32-bubble component. The first component
of size at least 16 appears at the saved 12 ms frame; all 32 belong to one component
at 50 ms. The matched Secondary OFF cloud has exactly zero displacement and no
topology events. No artificial motion, adhesion, fusion or cluster locking was added.

PBD corrections are symmetric, separating and preserve the centroid. For an
overlapping pair, their contribution to the change in global squared radius of
gyration is nonnegative, including the nonnegative squared-correction term; they do
not supply attraction. Cumulative correction/path ratios are nevertheless large in
packed clusters: 0.819 over Scene B and 1.049 over the 1 s run. These reflect repeated
reaction against acoustic compression, not a small correction to a resolved contact
hydrodynamic solution. Contact treatment can affect final packing and remains an
important calibration limitation. OFF stationarity and timestep statistics are
evidence about this implementation, not proof of physical accuracy.

## 3. Real PZT morphodynamics

Scene C uses the existing 2.5 MHz k-Wave complex basis, not an analytic substitute.
It begins with 32 dispersed bubbles in a 30 µm cloud near (0, 80, 40) µm. At 80 ms,
source amplitudes and phases change and Z− switches off. The run continues to 200 ms.
Source amplitudes are normalized to 850 Pa at the cloud center in each segment;
actual queried pressure spans 750.922–953.371 Pa. This normalization changes source
amplitudes without changing frequency or the stored field propagation solution.

At the identical switch state, the total force changes by 9.674% in Euclidean norm.
A separate unchanged-control continuation starts from exactly that saved state.
At 200 ms, the controlled trajectory differs from this continuation by 0.586 µm in
centroid, 8.523% in the shape-tensor norm and 5.747% in eigenvalue norm. Both trajectories
pass the fixed numerical-quality criteria. This demonstrates a computed control
response within the approximate model, not an experimentally calibrated effect size.

The controlled cloud contains a largest component of 18 bubbles at the switch and
all 32 by 90.8 ms in the saved frames. Accepted-step graph reconstruction records
36 merges and 5 splits; splits occur between 87.94 and 89.40 ms. These are changes
under the gap ≤ 0.1 R spatial threshold, sometimes transient, not bubble fusion or
physical breakup. The counterfactual is used only in plots/analysis; it is not
stitched into the continuous movie. Every saved partition agrees with reconstruction
of the per-substep edge log. Per-component tensors, membership and principal axes
are retained; near-degenerate axes have an explicit identifiability flag.

## 4. Model and valid domain

The archived adaptive APECSS table has eight valid leaves, with radius
1.779–1.781 µm, frequency 2.499–2.501 MHz and pressure 0.1–40000 Pa. F1 restricts
the actual trilinear cells to R = 1.78 µm and f = 2.5 MHz, without resampling or
extrapolation. The stored acoustic basis itself supports only 2.5 MHz. Zero drive
has exactly zero oscillatory response; a nonzero value outside the table is rejected.
Pressure and each volume harmonic retain the positive-exponential phasor convention;
harmonic h rotates by h times the pressure phase. Primary is the cycle average
−Re(V1 ∇P*)/2; Secondary sums separate h²-weighted harmonic correlations.

`smooth_near_field_exploratory` preserves the signed correlation C and uses
F = C Δx / (d² + a²)^(3/2), with a = (Ri+Rj)(1−q)² and
q = clip[d / (2(Ri+Rj)), 0, 1]. It is continuous, matches the far-field expression
above 2(Ri+Rj), and attenuates rather than sign-clips the near force. At contact
its magnitude is approximately 0.913 times the far-field value. It has no separate
velocity cap. `classical_far_field` is separately named and guarded below
5(Ri+Rj); this guard is a declared conservative numerical scope, not a validation
boundary established by experiment.

Translation uses Stokes drag with viscosity 0.001 Pa·s. The radial table retains
its own effective viscosity 0.002 Pa·s. Two relaxed block-Jacobi passes approximate
the reused Stage 4L local normal resistance, activated at h/R < 0.1 and regularized
at h/R = 0.001. A common velocity residual per current lubrication component enforces
sum(ζU) = sum(F), preserving all relative velocities within each component. It does
not lock a component or recenter the whole cloud. PBD uses diagonal Stokes mobility
weights (1/2 for the monodisperse bubbles), not a complete coupled hydrodynamic metric.

The initial prototype omitted that component force-balance condition and showed
2.85–3.50 µm of spurious uniform-field centroid drift. Its records are preserved and
explicitly retired in `PROTOTYPE_RETIREMENT.json`. The final Scene B drift is
1.19 × 10⁻²⁰ m. All final capabilities and timing claims below use the corrected
`fast_pzt_f1_v2_balanced_mobility` model.

## 5. Timestep selection and numerical quality

One representative formation scene screened 1/2/5/10/20/50 µs with four contact
passes. Larger candidates failed the predeclared 0.01 R residual-overlap limit.
Two targeted contact revisions used eight then twelve passes. After the force-balance
correctness fix, only 1 and 20 µs were rerun, both with twelve passes and unchanged
quality thresholds. The selected 20 µs case has Rg relative RMS difference 0.2953%,
mean-nearest-neighbor RMS difference 0.0616%, the same final component sizes, and the
same sampled formation time as its corrected 1 µs reference. Adjacency changes number
550 versus 498; these need not be identical for statistically equivalent packing.

The selected macro step has **one actual physical substep**, two mobility passes
and twelve local contact iterations. A 50 ms run has 2,500 macro / 2,500 actual
substeps; the 1 s run has 50,000 / 50,000. Contact iterations are disclosed work,
not concealed physical substeps. The maximum residual overlap is 13.082 nm in Scene B,
8.712 nm in Scene C and 16.880 nm across the batch, all below the original 17.8 nm
limit. The maximum allowed single contact-pass correction is 0.1 R = 178 nm.
All selected runs are finite, complete and in domain. Failed preliminary candidates
remain excluded; no quality threshold was loosened to obtain a pass.

## 6–7. Measured performance and long run

| Final-model case | Simulation wall | Including metrics + HDF5 | Process body including setup |
| --- | ---: | ---: | ---: |
| 200 bubbles, 5 ms | 0.352 s | 0.410 s | 0.762 s |
| 32 bubbles, 50 ms | 1.021 s | 1.067 s | 1.084 s |
| 32 bubbles, 1 s (actually completed) | 21.748 s | 22.027 s | 22.397 s |
| Real-basis 32 bubbles, 200 ms | 4.119 s | 4.565 s | 5.088 s |

Simulation wall time includes GPU execution, sampling synchronization and control
recombination; initial compilation and initial graph capture are excluded and
reported separately. Process-body time includes warm-cache initialization, dataset
metrics and HDF5, but excludes Python interpreter startup and final summary-JSON
writing. The field named `sampled_download_s` also includes waiting for queued GPU
work; it must not be interpreted as pure memory-transfer cost.

For N = 200, warm CUDA-event measurements are 13.496 µs for field/response,
12.800 µs for full pair force, 83.156 µs for one contact pass, and 1456.835 µs for
a whole step. They are separate benchmark graphs, not an exactly additive profile.
Contact is the leading measured cost. Reducing particle block size from 256 to 32
halved the preliminary step cost without changing the small-state formulas; this
layout was rechecked in the final model. Raw Python launch of 100 full steps takes
78.146 ms to enqueue and 141.604 ms including completion; empty synchronization costs
1.584 µs. Final graph enqueue for 100 steps is 5.341 ms.

The 200-bubble / 5 ms engineering target of 10 s is met. The completed **32-bubble**
1 s run meets the 60 s target. N = 200 over 1 s was not run: a naive short-run rate
extrapolation is about 70 s and is not a completed result. Performance depends on N,
packing, field choice, sampling and hardware; this is not a universal real-time claim.

## 8. Batch throughput

Eight independent seeds, each with 32 bubbles over 50 ms, were run serially and
together with exactly matching initial positions and sample times. All eight passed
numerical quality in both modes. Including sampled metrics and HDF5, serial time
was 8.591 s and batch time 1.363 s: **6.304×** more valid samples per measured GPU
wall hour, approximately 3,352 versus 21,133. This is the rate for 50 ms uniform-field
samples, excluding one-time startup, rendering and archival transfer. Including
per-run setup gives 8.711 s versus 1.382 s.

Atomic accumulation order is not bitwise deterministic. Across these seeds, the
largest same-ID difference was 2.638 µm after aggregation, while final largest
component size matched in all cases and the largest Rg relative difference was
0.841%. This is acceptable for the stated exploratory statistical task, not a
guarantee of reproducible microscopic packing or R3-level trajectory accuracy.

## 9. Use and remaining calibration

The qualified final datasets support exploratory aggregation/control studies and
small-batch software development within the stated domain. Near-field forces,
collective radial/acoustic scattering, packed-contact hydrodynamics and physical
parameter calibration still require work before quantitative prediction. The source
basis is an existing numerical field, not measured hardware calibration. The inherited
response table's validity flags and 5% production tolerance do not remove these limits.
No wall, background flow, added mass, history force or polydisperse batch is claimed.

Install/run commands are in [`scripts/fast_pzt/README.md`](../scripts/fast_pzt/README.md).
Four English figures and their source data are under `results/fast_pzt/f1/figures`.
The actual-radius, fixed-camera movie is under `results/fast_pzt/f1/visualization`:
1,001 actual frames at 60 fps, 0.2 ms physical time per frame, spanning 0–200 ms.
No interpolation, group enlargement or synthetic trajectory concatenation is used.
Large files are archived with local/server paths and SHA-256 in `ARTIFACT_MANIFEST.json`.

References used for implementation/scoping: [Warp installation](https://nvidia.github.io/warp/latest/user_guide/installation.html),
[Warp releases](https://github.com/NVIDIA/warp/releases),
[XPBD](https://mmacklin.com/xpbd.pdf),
[Doinikov & Bouakaz (2015)](https://doi.org/10.1103/PhysRevE.92.043001),
[Fan et al. (2020)](https://doi.org/10.1016/j.ultsonch.2019.104790).
The regularizer is not an implementation or reproduction of either near-field paper.
