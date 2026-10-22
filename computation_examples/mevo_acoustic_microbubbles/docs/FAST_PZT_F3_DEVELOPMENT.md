# F3 — PZT Shape Controllability & Numerical Reproducibility

Date: 2026-10-08. Branch: `dev/fast-pzt-aggregation-f3-20261008`.
Base F2 archive: `6655867b772484a82ea9a7c2afd491da02bdd378`.
Results: `results/fast_pzt/f3/`. Physical scope: **EXPLORATORY**.

## Scientific outcome

The fixed-order GPU implementation reproduces the three tested 120 ms trajectories
bit for bit on the same device. PZT phase and amplitude/enable changes produce
resolved centroid motion. They also produce brief shape-eigenvalue responses, but
neither control establishes a sustained, robust shape response under the frozen
criteria and separately measured sensitivity envelopes. No persistent split or merge
occurs. Determinism has improved the causal comparison; it has not validated the
physical accuracy or perturbation robustness of dense-cluster rearrangement.

| Capability | F3 outcome |
|---|---|
| PZT differential driving | `FIELD_DIFFERENTIAL_DRIVE_LIMITED` for direct Primary drive; measurable Total-force modulation exists |
| Numerical reproducibility | `PASS_SAME_DEVICE_BITWISE`, three repeats; no cross-device promise |
| Translation control | Resolved: phase 0.2500 µm, amplitude/enable 1.3390 µm relative to unchanged at 120 ms |
| Shape control | `TRANSIENT_SHAPE_RESPONSE`; both `NO_ROBUST_ENDPOINT_SHAPE_RESPONSE` |
| Connectivity | `NO_PERSISTENT_SPLIT`; no component merge or member exchange |
| Numerical quality | All 12 trajectories pass the inherited checks; this is not physical validation |
| Performance | N=200, 5 ms: 0.7074 s atomic; 0.7270 s stable |
| Remaining sensitivity | `MORPHOLOGY_NUMERICALLY_SENSITIVE` |
| Physical scope | `EXPLORATORY`; no new calibration |

`FIELD_CONTROL_LIMITED` describes the tested direct spatial drive and the absence
of established sustained shape control. It does not mean that changing PZT controls
has no force effect. In particular, the selected amplitude control changes the
existing Secondary-force pattern appreciably.

## Frozen model and experimental order

The archived F2 mature state at absolute time 102.2 ms supplies all 32 actual IDs and
positions. Every unperturbed continuation starts from exactly those bytes; the
initial cloud is not regenerated. The saved mature-state file hash is
`eb1ba215ee79bf8d8175b91f99360f72524ade74fe0e60e3c26872d7d2f1a8f4`.

F2 pair-feedback near field, F1 central PBD contact, float64 Warp full-pair force,
two local mobility iterations, twelve contact iterations, 10 µs macro/physical step,
one actual substep, R=1.78 µm, f=2.5 MHz, eight harmonics, response table, and four
PZT basis fields are unchanged. CUDA Graphs and independent particle degrees of
freedom remain. There is no vessel, imposed flow, adhesion, coalescence, speed cap,
or cluster locking in these F3 experiments. The rejected F2 tensor contact is not
used. No R3 high-fidelity job or old stopped task was resumed.

The committed plan and `configs/fast_pzt/f3/protocol.json` preceded the trajectories.
Six fixed-state candidate slots were examined first; reduction checks and three
repeats per mode followed. Only the selected phase and amplitude controls entered
the dynamics comparison. Two predeclared initial perturbations and one N200/5ms
run per mode complete the bounded campaign. There was no search for attractive
final shapes and no post-trajectory control tuning.

## A. Actual field and differential forces

All force components, complex pressures, gradients, response harmonics, velocities,
component roots, and control-minus-unchanged centered differences are retained in
`driving_v3/FORCE_FIELD_DATA.h5`. `CONTROL_FORCE_DIAGNOSTICS.json` adds the tensor,
fit, grid-consistency, Secondary-scaling and angular-balance diagnostics.

The reference force is fixed at Stokes drag for 1 µm/s:
F_scale = 6πµR × 1 µm/s = 3.35522e-14 N = 0.0335522 pN.
DP is RMS of centered Primary force; DS is RMS Secondary force. The ratio
D=DP/max(DS,F_scale) is descriptive, not a calibrated deformation threshold.

| Control | DP (pN) | DS (pN) | D | RMS centered control-induced Primary (pN) | RMS centered control-induced Total (pN) |
|---|---:|---:|---:|---:|---:|
| Unchanged | 0.008056 | 165.8202 | 4.8580e-5 | 0 | 0 |
| Phase only | 0.008238 | 165.7485 | 4.9701e-5 | 0.002128 | 0.694190 |
| Amplitude/enable | 0.005819 | 159.8062 | 3.6415e-5 | 0.004407 | 7.129805 |

The unchanged mean Primary vector is approximately
[-4.7337e-17, -1.70287e-13, -8.73498e-14] N. Its spatial nonuniformity is 4.21%.
The selected controls' differential Primary increments are only 1.28e-5 and
2.66e-5 of baseline DS. Direct differential Primary actuation is therefore small
in this measured configuration, not merely inferred from wavelength or grid size.

Total-force increments are 0.419% and 4.300% of baseline DS. The amplitude candidate
passes the predeclared **instantaneous** strong-drive diagnostic through Total force
and the mobility strain score. A fixed-position least-squares decomposition
S_control = α S_unchanged + residual gives α=0.963466: 72.25% of the amplitude
control's Secondary-change energy is global scaling of the existing, spatially
nonuniform Secondary pattern. Its residual RMS is 3.7544 pN. This provides a real
route to modifying compression and rearrangement, but does not establish selective
shape actuation. For phase control the global-scaling fraction is only 1.11%, with
a much smaller 0.69035 pN residual.

### Spatial resolution and interpolation

The 32 centers span 9.53 × 10.81 × 12.05 µm; including actual radii the extent is
13.09 × 14.37 × 15.61 µm. They occupy only cells [101,18,18] and [102,18,18]
(14 and 18 centers) on the 100 µm grid. Both actual particle positions and a
5×5×5 ROI padded by 2R are evaluated.

Production gradients are nodal `numpy.gradient` values interpolated separately
from pressure. They are not the derivative of trilinearly interpolated pressure.
Centered pressure finite differences remain inside each particle's cell, use
at most 0.05 µm offsets, and agree after halving the offset at approximately 1e-12
relative error. Nevertheless, stored-gradient versus pressure-derivative L2
differences are **24.2–38.0%** across eligible candidates. The corresponding
Primary-force discrepancy is also substantial. No production gradient was changed
in F3. This is an unresolved field-discretization inconsistency, not an estimated
experimental error bar and not proof that either interpolation is physically exact.

### Kinematics before contact

The position design matrix has rank 3 and condition number 1.2674. Free-Stokes,
two-iteration local mobility, and mean-force-balanced velocities have separate
L/E/W fits; degenerate fits withhold tensors. Baseline centered velocity RMS falls
from 4.9423 mm/s in free Stokes to 1.8516 mm/s after local mobility. These are
pre-contact predictions; contact strongly alters the accepted motion.

The centered control-induced mobility velocity is 3.8823 µm/s for phase and
74.1069 µm/s for amplitude. Their deviatoric symmetric-gradient scores are
0.25519 s⁻¹ and 2.05694 s⁻¹. The relative affine-fit residuals are **95.99%** and
**56.23%**: the score summarizes a fitted component and does not describe all local
motion or a validated strain rate. Full L, E, W and residuals remain in JSON.

An additional instantaneous audit compares acoustic torque to ζΣr×U before contact.
At baseline the torque residual has norm approximately 6.66e-16 N m, while the
acoustic torque is approximately 1.19e-19 N m. Two truncated local mobility passes
plus a common mean-force correction do not enforce the angular balance of an
exact symmetric normal-resistance solve. This independent approximation limits
rotation interpretation. It does not identify the sole cause of later rearrangement,
and was not corrected by changing the frozen physics.

## B. GPU reduction and reproducibility

`f3_engine.py` subclasses the existing F2 engine. `fast_atomic` delegates the
original force/component implementation. `stable_reduction` writes each unique
pair contribution once into a persistent GPU buffer, then sums partners in
ascending particle-ID order. Component corrections are summed in ascending
member-ID order. The pair expression, precision, contact/mobility kernels,
integration clock and integer component graph are unchanged. Diagnostic
`stable_pairs` and `stable_components` modes isolate the two reductions.

Eight fixed-state evaluations per mode show up to 7.7548e-26 N variation in atomic
Secondary/Total force while Primary, pressure and harmonics remain identical.
Stabilizing components alone retains that variation. Stabilizing pair accumulation
removes fixed-state variation in this sample, including downstream component
corrections. Component atomics alone did **not** vary in these eight samples; their
independent contribution is not experimentally established. This evidence supports
pair accumulation order as an observed seed of trajectory divergence, without
claiming that all dense-cluster sensitivity originates in atomics.

The fixed-state NumPy F2 pair formula agrees within rtol=1e-11, atol=1e-24 N.
Canonical component membership agrees with a CPU graph. Stable CUDA Graph and
direct launches give identical positions and sorted edge events for 50 steps.
The dispersed OFF test remains stationary and an off-grid 17 µs switch occurs
at the exact registered time.

| Three identical-input runs | Atomic | Stable reduction |
|---|---:|---:|
| Maximum repeat-pair Rg difference, any time | 1.1685% | 0 |
| Maximum repeat-pair eigenvalue difference, any time | 8.2044% | 0 |
| Maximum fixed-ID aligned internal RMSD | 1.90966 µm | raw states identical |
| Maximum disagreeing neighbor edges | 65 | 0 |
| Simulation wall times (s) | 5.4887 / 5.4667 / 5.4615 | 5.6018 / 5.6042 / 5.6022 |

All saved stable positions and sorted raw edge events agree bit for bit. SVD-based
alignment of a trajectory with itself leaves about 1e-20 m arithmetic residual;
that is not a trajectory difference. The envelope is the maximum of all three
repeat-pair comparisons, not a statistical confidence interval. These measurements
do not replace F2's historical two-run sensitivity values or imply cross-GPU,
cross-compiler or cross-block-layout bitwise portability.

Two separate initial-condition experiments add a fixed, zero-mean displacement
direction of RMS 1 pm and 1 nm. Both initial states pass overlap checks without
projection. Their endpoint internal RMSDs are 3.42423 and 3.43018 µm relative to
the unperturbed stable trajectory. Endpoint Rg differences are 2.7421% and 1.2620%;
eigenvalue differences are 14.6951% and 4.0841%. Identical-input reproducibility
therefore coexists with large sensitivity to actual initial perturbations.

## C. Finite control selection

Element order is Y+, Y−, Z+, Z−. Unchanged amplitudes are all 151.340514 Pa,
with all elements ON and all phases zero. Selected phase control preserves those
amplitudes/enables and uses phases [30°,45.67023°,0°,0°]. Selected amplitude control
uses [0,246.257925,0,246.257925] Pa, phases zero, only Y− and Z− ON.

The six slots are unchanged, the F2 phase control, one prescribed Z-pair analytic
phase candidate, and the three existing F2 amplitude patterns. The prescribed
Z+=30° ROI-RMS matching equation has no real Z− root. That slot is recorded as
ineligible; its Z−=0 field diagnostic does not become an extra dynamic candidate.
The eligible amplitude scores are 0.99968, 1.21018 and 2.05694 s⁻¹; the third wins.
The remaining eligible phase candidate is selected. This bounded comparison is
not a global optimization claim.

Phase particle/ROI RMS-pressure changes are 0.01938% and 0.000650%; amplitude
changes are 1.8052% and 1.6793%. All are retained explicitly; the amplitude case
is not exactly isobaric. Baseline particle RMS is 860.243 Pa and ROI RMS 859.696 Pa.
Eligibility uses ±10% RMS tolerance at both sets of points, response/basis domain,
and the conservative F2 campaign source-amplitude bound 291.237688 Pa. This bound
is **not** a measured hardware rating. Matching F2's 27 diagnostic points alone
would not guarantee equality in this actual small ROI.

## D. Actual shape response

The unchanged stable repeat supplies the matched reference. Each selected control
continues for 120 ms, ending at absolute 222.2 ms, with 601 saved states at 0.2 ms
spacing. No ID reassignment, best-looking camera selection or radius adjustment
is used. Fixed-ID metrics and instantaneous component metrics are stored separately.

| Endpoint versus unchanged | Phase only | Amplitude/enable |
|---|---:|---:|
| Centroid distance | 0.250023 µm | 1.339023 µm |
| Relative Rg difference | 0.03654% | 1.59355% |
| Relative eigenvalue-vector difference | 0.15994% | 3.92346% |
| Fixed-ID rigid-aligned internal RMSD | 2.84122 µm | 2.70764 µm |
| Disagreeing neighbor edges | 99 | 100 |
| Final aspect ratio | 1.22319 | 1.23352 |

The unchanged final aspect ratio is 1.22087. Kabsch angles of 168.1° and 125.7°
coexist with substantial internal rearrangement and do not establish controllable
rigid rotation. Principal-axis orientations are masked when the relative eigengap
is below 0.1; changing a nearly degenerate axis is not evidence of rotation.

Shape detection requires the frozen absolute scales (2% Rg, 5% eigenvalue-vector,
0.1R aligned RMSD), greater than 3× the stable repeat envelope, and comparison
with the pointwise maximum of historical F2 same-mature-state near-field/mobility
contrasts and the two F3 initial-perturbation contrasts. Historical variants were
not rerun. This is a conservative **descriptive sensitivity envelope**, not a
calibrated model-error bound or confidence interval.

Phase eigenvalue contrast reaches 7.1439%; its longest interval above all criteria
is 15.0–17.6 ms after switching (2.6 ms). Amplitude reaches 13.7164%; its longest
such interval is 70.8–79.6 ms (8.8 ms). Its Rg response exceeds the criteria for
only 0.8 ms. Neither reaches the registered 10 ms sustained duration. Internal
RMSD does not exceed the combined sensitivity envelope. Endpoint contrasts are
below that envelope (3.3905% Rg, 14.9217% eigenvalues, 4.16896 µm internal RMSD).

Both controls are therefore `TRANSIENT_SHAPE_RESPONSE` and
`NO_ROBUST_ENDPOINT_SHAPE_RESPONSE`. The comparison detects transient numerical
model responses, not validated lasting physical elongation or rotation. All 32
members remain in one component. Raw pair-edge flips (unchanged 670, phase 633,
amplitude 942) are neighbor rearrangements, not component splits or physical fusion.
There are zero persistent component splits/merges and no membership exchange.

## E. Quality and performance

All six repeats, two control branches, two initial perturbations and two performance
runs finish with finite positions, zero domain/error bits, exact intended clocks,
and the original quality bounds: accepted overlap ≤17.8 nm and per-pass contact
correction ≤178 nm. The mature reference/phase/amplitude maxima are respectively
10.676/11.052/9.798 nm accepted overlap and 18.005/16.300/16.021 nm per-pass correction.
The 1 pm perturbation has the largest accepted overlap among these branches,
12.701 nm. No threshold was loosened.

These checks alone do not establish precise dense-cluster dynamics. Summed contact
correction magnitude divided by summed mobility-predicted path is approximately
1.045, 1.047 and 1.037 for reference/phase/amplitude. Relative to the free-Stokes
path it is 0.436, 0.443 and 0.411. These sums are not net physical displacements;
their large magnitude shows that dense contacts strongly influence accepted
motion. That limitation, the angular-balance residual, initial sensitivity and
field-gradient inconsistency remain visible beside the numerical-quality PASS.

| N=200, 5 ms, 500 physical steps | Fast atomic | Stable reduction |
|---|---:|---:|
| Simulation loop wall time | 0.70736 s | 0.72701 s |
| Loop + metrics + HDF5 (`dataset_wall_s`) | 0.76255 s | 0.78195 s |
| Pair-force stage, CUDA-event mean | 22.2925 µs | 42.6598 µs |
| Whole step, CUDA-event mean | 1514.53 µs | 1540.64 µs |

Stable/atomic loop-wall ratio is 1.02778 (+2.78%). The earlier F2 loop measurement
was 0.72433 s. The pair stage costs more because it writes a full persistent pair
buffer and performs an additional ordered reduction. Its near doubling does not
double total cost: contact/mobility and other work dominate. CUDA-event measurements
use 100 warmed graph repetitions and are not additive component decompositions.
`sampled_download_s` includes waiting for queued GPU work; it must not be interpreted
as pure transfer time. These are single matched runs, not a statistical benchmark
or new batch-throughput validation. Stable reduction is the F3 causal-study default;
fast atomic remains available for high-throughput exploratory generation. The extra
pair buffer scales as 24 × batch × N² bytes and is not a large-N scalability proof.

## Figures, video and identity

`figures/Fig_F3_1_differential_driving.{png,pdf}` shows actual spheres, Primary and
centered control-induced Primary vectors, force scales and pre-contact relative
velocities. Arrow scales (µm per pN) are explicitly labeled and stored in
`VECTOR_SCALES.json`; spheres always retain the actual radius. Surface-anchored
arrows improve visibility without shifting particle centers.

`Fig_F3_2_numerical_reproducibility` compares three repeats per mode.
`Fig_F3_3_PZT_control_response` separates translation, transient intervals,
endpoint differences and sensitivity envelopes. Single-sample window markers have
a 0.2 ms display width; their inferred observed duration remains zero in JSON.
`Fig_F3_4_actual_morphology` uses identical cameras and coordinate bounds for the
initial state and three endpoints. PNG/PDF and plotted arrays are retained.

`visualization/F3_mature_PZT_comparison_actual_radius_60fps.mp4` is a diagnostic
comparison, not a claimed shape transformation. It contains all 601 actual saved
frames, 1800×800 pixels, 60 fps, 10.0167 s encoded duration, with true 0.2 ms
physical increments and absolute times 102.2–222.2 ms. There is no interpolation,
particle subset, cluster magnification or moving camera. Decoder validation,
time mapping and input paths are in `VIDEO_VALIDATION.json`.

Scientific executions and archive commits are distinct:

| Successful F3 job | Scientific execution commit |
|---|---|
| checks_v2 | `df7637fab8b2e59cb2012e364b890a524dcc9b92` |
| driving_v3 | `3cc83ec0bf70c15569521772f232b4d5a0f865e6` |
| repeats_v1 | `f7c3469d634560d539cda6ca2379f130059e8022` |
| controls_v1, performance_v1 | `58ea583a78a09938b8dc6ddcf773d9df64cb37ff` |

F2's archive is `6655867...`; its historical science remains attributed to
`31e66ce...` (core), `cd7d6f5...` (sensitivity), `bfd228a...` (mature),
`8feabe0...` (performance), and `03655ec...` (phase). Per-run identities retain
full hashes. Offline renderer file hashes and the later delivery-source commit
are separately indexed; none is substituted for a trajectory's actual execution.

Failed F3 startup attempts (`driving_v1`, `checks_v1`) contain no simulation;
an import re-export was fixed. `driving_v2` stopped when the prescribed analytic
phase equation had no feasible root; the candidate was then explicitly marked
ineligible. Logs remain archived. Initially misplaced F3 startup logs were moved
from the F2 root; no existing F2 artifact was overwritten. Server checks confirm
all 243 F1 and 174 F2 manifest artifacts retain their hashes.

The response table and basis SHA-256 values are respectively
`5d53cb2ef42fa64591f45e3e7de5db106d6cbe74794c8a9abfa88a3354bec983` and
`39c777eb3476fbc3878b48f2adac2b9dc3a81be1090f40756fc80ff3980e0266`.
The 33 MB basis is indexed at its preserved local/server locations rather than
duplicated in Git. All new raw trajectories and the short MP4 are below the 4 MiB
per-file archive limit. `ARTIFACT_MANIFEST.json` records sizes and hashes.
See `scripts/fast_pzt_f3/README.md` for reproduction commands and environment use.

The supported next scientific questions are field/gradient consistency and robust
dense-contact mobility/shape inference. F3 does not silently repair either or
claim experimental quantitative calibration, stable controllable rotation, lasting
elongation, or physical split/merge.
