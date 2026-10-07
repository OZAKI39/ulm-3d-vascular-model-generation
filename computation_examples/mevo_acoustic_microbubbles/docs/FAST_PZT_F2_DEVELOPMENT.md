# F2 — Near-Field Interaction & PZT-Controlled Cluster Morphodynamics

F2 development and the bounded GPU campaign are complete. The selected backend is
**F2 pair-feedback near field + F1 central PBD, 10 µs**. It preserves Warp graphs,
persistent arrays, eight harmonics, full pairs, independent bubble coordinates and
batching. The new anisotropic contact candidate is **rejected as the default**.
PZT changes move a mature cluster and produce descriptive transient shape responses;
lasting endpoint shape control is **not established**. Physics remains **EXPLORATORY**.
Read the separate gates in [`FINAL_GATE.json`](../results/fast_pzt/f2/FINAL_GATE.json).

## 1. Identity, preservation and actual execution

Branch: `dev/fast-pzt-aggregation-f2-20261008`, based exactly on F1 archive
`6de76d508ea088e137069abf058b86b7f6d47a1c`. The final archive commit is the branch
commit containing this report; it must not replace historical execution identities.
Historical F1 science remains `9829156893339efff713483f3125848b84adbc98`, with matched
controls at `3721bebd597164262691af766f723d245e043116`.

| F2 evidence | Actual execution commit |
| --- | --- |
| Core formula/GPU checks and two-bubble trajectories | `31e66ce73c26e421101c335009f20a9462179a40` |
| Ordered force/contact/step comparison | `cd7d6f5ca4c8a3298b5e7aedf2f64642f9c85316` |
| Contact disposition and mature-state controls | `bfd228ab69fe7fe7a3ba29c43b06cf53d0f4f3c2` |
| 200-bubble performance regression | `8feabe08c664a190cdea719dbd43aee96e3a4856` |
| RMS-matched phase-only continuation | `03655ec93f2cb518a428d9db496a5312258d1575` |

Each run also stores source/input SHA-256 in `IDENTITY.json`. All six original F1
simulator modules are byte-identical to the exact F1 archive; see `F1_PRESERVATION.json`.
Stage 4/4L/4M, Stage 5A-R3, F1 results and old failures are preserved. No old long
job was resumed. New results reside in the isolated local worktree
`/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f2_20261008/computation_examples/mevo_acoustic_microbubbles/results/fast_pzt/f2`
and server `/workspace/mevo_fast_pzt_f2/results/fast_pzt/f2`.

The F1 baseline was inspected, not rerun as a full campaign: 20 µs, one substep,
two mobility iterations, twelve contacts, a complete 32-bubble component by 50 ms,
200 bubbles / 5 ms in 0.351664 s, 32 bubbles / 1 s in about 21.748 s, and historical
8-sample batching speedup about 6.304×. Its correction/path ratios about 0.819 and
1.049 were recognized as a packing-sensitivity concern, not evidence of attraction
created by contact.

## 2. Core changes and physical hypothesis

New modules are `nearfield.py`, `contact_models.py`, `f2_kernels.py`, `f2_engine.py`
and `morphodynamics.py` under `src/mevo_mb/fast_sim`. F2 subclasses the existing
engine and overrides pair/contact stages. The F1 backend remains independently
available; choosing both F1 model names delegates to its original kernels.

For the positive-exponential phasor convention, monopole scattered pressure is
`p_sc,j = -K V_j`, with `K = rho omega² / (4 pi d)`. Let the local table secant
compliance be `chi_i = V_i1^0 / P_i`. Solve the two complex scalar equations:

```
V_i + K chi_i V_j = V_i1^0
V_j + K chi_j V_i = V_j1^0
V_i = (V_i1^0 - K chi_i V_j1^0) / (1 - K² chi_i chi_j)
```

Only the signed fundamental correlation is replaced:
`C1 = (1-w) Re(V_i1^0 V_j1^0*) + w Re(V_i V_j*)`.
The cubic smoothstep is one at `d <= Ri+Rj`, zero at `d >= 4(Ri+Rj)`; it is continuous
with zero endpoint derivative. Harmonics 2–8 keep their original signed, h²-weighted
contributions. Both candidates retain F1's geometric soft-core regularizer, which
becomes classical beyond `2(Ri+Rj)`. The complete F2 interaction becomes the classical
uncoupled far-field expression beyond `4(Ri+Rj)`.

There is no fitted attraction multiplier or sign clipping. General pair correlations
can change sign through phase feedback; identical in-phase and counterphase formula
checks retain attraction and repulsion. At 1 kPa and geometric contact, the identical
in-phase response enhancement is about 1.60 in fundamental amplitude (about 2.56 in
fundamental force correlation), arising from this declared closure, not a cluster-size
fit. The formula is **not calibrated**.

The model is a force-only pair approximation. Different pairs imply different effective
radial responses for the same bubble; there is no unique collective response. Primary
force remains the original uncoupled table Primary. Higher-harmonic scattering, nonlinear
radial self-consistency and multipole translation are absent. Zero-drive compliance
uses the table's lowest valid-pressure secant, while all-OFF response remains zero.
Closure singularities and out-of-domain coupled effective pressures are explicit errors.

[Doinikov & Bouakaz (2015)](https://doi.org/10.1103/PhysRevE.92.043001) motivates the
limits of far-field theory; its publisher abstract was inspected, not its full
unavailable PDF. [Fan et al. (2020)](https://doi.org/10.1016/j.ultsonch.2019.104790)
is a radial-model scope reference, not an implemented aggregation force. Neither
paper is claimed as reproduced. [XPBD](https://mmacklin.com/xpbd.pdf) informed the
constraint-metric interpretation; inertial mass weights were not transplanted into
this overdamped model.

## 3. Near-field sensitivity: actual matched experiments

Two GPU trajectories began at 8 µm separation with identical 1 kPa input. Formula
checks compared Warp to an independent complex 2×2 NumPy solve at 82 distances,
including artificial opposite-phase responses, continuity and reciprocity. Both
trajectories approached contact without violating quality. The F2 pair approaches
faster; signed curves, sampled relative speed and raw pair trajectories are retained.

The 32-bubble comparison reuses Scene B seed 20261007, 20 µm cloud, uniform 1 kPa,
50 ms duration, the same positions/contact and **matched 10 µs**:

| Metric | F1 near field | F2 near field |
| --- | ---: | ---: |
| First sampled Nmax ≥ 16 | 12.0 ms | 10.5 ms |
| First sampled Nmax ≥ 29 | 14.5 ms | 13.5 ms |
| First complete 32-bubble component | 16.0 ms | 14.5 ms |
| Final connected components / largest size | 1 / 32 | 1 / 32 |
| Maximum accepted overlap | 7.206 nm | 14.102 nm |
| Sum of contact-correction lengths | 1.519 mm | 3.070 mm |
| Correction / mobility-predicted path | 0.808 | 0.886 |
| Simulation wall | 2.046 s | 2.143 s |

Final Rg differs by **0.562%**, shape eigenvalues by **5.504%**, while fixed-ID,
rigid-aligned RMSD is **3.824 µm**. Raw edge changes are 514 versus 900. The half-group
formation time improves 12.5%, while full-component first appearance improves 9.375%.
`formation_time_s` in run summaries is explicitly the inherited first Nmax ≥ 16
metric, not completion or maturity. Derived first-hit fractions disambiguate it.

Coarse Rg/eigenvalue/formation differences fall below the registered 5%/10%/20%
sensitivity flags, but individual packing is clearly not invariant. This single-seed
comparison does not identify which packing is physically correct. Parameter values
were not tuned to force agreement.

## 4. Contact and mobility: sensitivity and rejected candidate

The candidate freezes `Mi = [I + sum cij nij nijᵀ]^-1` at the predicted state using
F1's normal resistance. With `Hij = (Mi+Mj)/2`, an overlapping pair receives
`delta_i = 0.5 overlap Hij n / (nᵀHij n)`, `delta_j = -delta_i`, with 0.5 Jacobi
relaxation and twelve passes. It preserves the center of mass and has separating
normal work; one contact pass cannot reduce global squared Rg. It adds no adhesion,
friction, fusion, cluster freezing or QP.

However, equal-and-opposite anisotropic corrections are not central. An independent
fixed irregular four-bubble configuration gives angular correction moment
`[-3.09323e-16, ~0, ~0] m²` and an infinitesimal fitted rotation of
**3.11896e-5 rad per pass**. This is a nonzero algorithmic angular bias. The candidate
is preserved for diagnostics but **rejected for mature-cluster causal control claims**.
Its initial numerical selection at 5 µs is superseded by `adoption_v1/SELECTED.json`;
no historical selection file was rewritten.

Matched contact comparison at **5 µs**, with F2 force fixed:

| Metric | F1 central PBD | Mobility candidate |
| --- | ---: | ---: |
| First Nmax ≥ 16 | 10.5 ms | 10.5 ms |
| First complete component | 15.0 ms | 14.5 ms |
| Final components / size | 1 / 32 | 1 / 32 |
| Maximum accepted overlap | 6.909 nm | 5.893 nm |
| Sum of correction lengths | 2.814 mm | 2.713 mm |
| Correction / mobility-predicted path | 0.875 | 0.843 |
| Correction / free-Stokes path | 0.334 | 0.328 |
| Simulation wall | 4.273 s | 8.665 s |

Final Rg differs **2.341%**, eigenvalues **12.452%**, tensor Frobenius norm **15.817%**,
and ID-aligned RMSD **3.669 µm**. Raw edge changes are 672 versus 960. The eigenvalue
flag exceeds the registered 10% criterion: **CONTACT_MODEL_SENSITIVE**. A lower
correction ratio alone does not qualify the candidate. The matched Secondary OFF
case remains exactly motionless, with 32 singleton components: contact does not
spontaneously attract the dispersed cloud.

A same-mature-state unchanged-control branch with four mobility iterations instead
of two gives a 3.390% final Rg difference, 7.186% eigenvalue difference and 3.779 µm
ID-aligned RMSD. Its correction/path ratio is 0.911 instead of 1.021, and wall time
7.180 s instead of 5.474 s. Thus coarse size is relatively stable but internal packing
and rearrangement remain sensitive. Linear force balance is not proof of local
hydrodynamic or angular accuracy. Neither central PBD nor the truncated local mobility
has been validated as a quantitative dense-cluster packing model.

## 5. Mature-cluster PZT experiment

The real saved 2.5 MHz basis forms a 32-bubble cloud (seed 20261008, radius 30 µm,
center [0,80,40] µm). Baseline source amplitudes are 151.3405 Pa for Y+, Y−, Z+, Z−,
all phases zero and all enabled. The first state meeting the frozen maturity rule is
**102.2 ms**, containing **all 32 IDs**. The rule requires ≥90% membership, unchanged
connected membership throughout 10 ms, Rg window range ≤2%, and the original overlap
budget. Raw accepted-step graph changes verify the window. Formation was bounded
at 160 ms and cost 6.763 s; only the actual prefix through the qualifying state is
used for the control branches and movie.

All six 120 ms continuations start bitwise at the saved mature coordinates. This is
a complete dynamical state for the memoryless overdamped integrator, not for an
inertial/history model. Unchanged, repeated unchanged, phase-only, amplitude/enable,
F1-near-field sensitivity and four-mobility-iteration sensitivity are retained.
All use 10 µs and complete 12,000 actual steps.

Control selection uses only field diagnostics at 27 fixed ROI points, before the
respective controlled trajectory. The amplitude choice is [0,246.2579,0,246.2579] Pa,
Y−/Z− enabled, all phases zero. Its ROI RMS pressure matches baseline 843.7695 Pa;
normalized complex-gradient contrast is 1.0294.

The initial three phase patterns changed ROI RMS by 16.7–27.5%, violating the
unchanged ±10% comparability criterion. They are preserved as ineligible diagnostics.
A documented analytic repair sets Y+ to +30° and solves the ROI RMS equation for Y−,
choosing the root with smallest absolute phase. The resulting phases are
**[30°,45.67°,0°,0°]**, with all source amplitudes/enables unchanged. ROI RMS matches
exactly; gradient contrast is 0.5179. No particle motion or force parameter was chosen
to manufacture a target shape. Matching regional RMS does not imply matching local
pressure everywhere; complex pressures and gradients at every diagnostic point are saved.

### Endpoint comparisons after 120 ms continuation

All values below compare with the same-state unchanged branch, not with the initial
cloud or an unrelated trajectory:

| Contrast | Amplitude/enable | Phase only | Same-control repeat |
| --- | ---: | ---: | ---: |
| Centroid distance | 1.339 µm | 0.250 µm | 0.0000173 µm |
| Relative Rg difference | 1.006% | 1.156% | 1.120% |
| Relative eigenvalue difference | 17.254% | 15.740% | 15.600% |
| Aligned internal RMSD | 3.995 µm | 4.005 µm | 3.203 µm |
| Kabsch best-fit rotation | 121.38° | 41.81° | 153.33° |

The centroid response is resolved above the repeat and tested near-field/mobility
comparisons. Its distance is a vector difference between trajectories: it is not
the difference of scalar travel distances shown in Figure F2-3c. Amplitude control
moves approximately in the opposite y direction to unchanged control, explaining
why the inter-trajectory distance exceeds either individual travel distance.

Endpoint size/eigenvalue/internal/rotation contrasts fail the registered requirement
of both an absolute effect threshold and more than three times the same-control
repeat difference. **NO_SIGNIFICANT_SHAPE_RESPONSE applies to the endpoint**, not to
all times. One repeat supplies a descriptive numerical envelope, not a confidence
interval. Atomic accumulation is not bitwise deterministic, and dense packing amplifies
small differences; the large repeat rearrangement prevents a reliable claim of
controlled microscopic rearrangement. Contact robustness remains unestablished.

### Transient responses, rotation and topology

Time-resolved screening retains the same thresholds. Eigenvalue contrast exceeds
both the repeat and tested near-field/mobility envelopes in short intervals: the
longest sampled span is **5.4 ms** for amplitude control (31.6–37.0 ms after switching)
and **5.6 ms** for phase control (36.0–41.6 ms). These are descriptive transient shape
responses. Fixed-ID aligned internal contrast never exceeds the tested model envelope.
No multiple-comparison significance or lasting elongation is inferred.

Final aspect ratios are 1.055 unchanged, 1.254 amplitude and 1.223 phase, versus 1.267
at maturity. Principal axes are masked when the leading relative eigenvalue gap is
below 0.1. Large Kabsch angles coexist with multi-micrometer residual deformation;
they cannot be interpreted as isolated rigid rotation or proof of PZT orientation
control. Both rotation matrices and aligned residuals remain available.

All mature branches keep the same 32-member connected component. **No member exchange,
raw component merge/split, or persistent component merge/split occurs.** Pair-neighbor
changes do occur: raw flips number 1479 unchanged, 979 amplitude, 633 phase. Final
neighbor Jaccard indices relative to maturity are 0.326, 0.346 and 0.508 respectively.
Among changed-pair episodes, linked episodes ending within 2 ms number 165, 108 and
93; brief unlinked episodes number 314, 203 and 131. These threshold crossings affect
internal edges, not component membership. They are never described as real bubble
fusion or breakup. Raw events and left/right censoring are preserved independently.

## 6. Numerical quality and step selection

The 0.01R accepted-overlap (17.8 nm) and 0.1R single-contact-pass correction (178 nm)
limits were unchanged. The nominal 20 µs F2-near/F1-contact case fails at **36.054 nm**.
The new contact candidate fails at **130.882 nm** (20 µs) and **58.317 nm** (10 µs).
Only these failures triggered the ordered 10/5 µs comparisons; there was no 1–50 µs
sweep, 1 µs fallback, threshold relaxation or velocity clipping.

Selected 10 µs F2 near-field / F1 contact passes. The 32-bubble uniform comparison
has maximum accepted overlap 14.102 nm. Mature formation and all six continuations
pass, with maximum accepted overlap up to 11.0525 nm. The maximum single-pass correction
across selected runs is recorded explicitly in `FINAL_GATE.json`; all are below 178 nm.
There are no nonfinite positions, response/field domain errors or closure errors.
PZT OFF and Secondary OFF give no spontaneous motion/aggregation. The noncommensurate
17 µs control-event check confirms the exact event clock. Six small GPU checks and
six offline invariance/data checks pass. They do not imply physical calibration.

For unchanged/amplitude/phase mature continuations, correction / mobility-predicted
path is **1.021 / 1.032 / 1.047**. Relative to free-Stokes path, it is
**0.399 / 0.405 / 0.443**. The two denominators are deliberately distinct: the former
preserves F1's metric, while the latter explicitly measures the force/zeta free path.
Contact remains a large part of dense dynamics; reducing its ratio was not an objective
that could override morphology or invariant failures.

## 7. One short performance regression

Same RTX 4090, same 200 initial positions/seed/cloud/1 kPa input as F1, 5 ms physical
horizon, selected F2 backend:

- Simulation wall **0.724332 s**, versus F1 **0.351664 s**: **2.0597×**.
- 500 macro steps and 500 actual substeps, versus F1's 250; per-step wall ratio is
  about **1.030×**. The main total-cost increase is the necessary halved time step.
- Dataset cost including metrics/HDF5 **0.779673 s**; process body **1.153896 s**.
- Accepted overlap **0.004206 nm**, no domain/finite errors; final largest component
  five bubbles. This short regression is not a completed 200-bubble aggregation claim.

CUDA-event measurements in warmed repeated-kernel graphs: field/response 13.599 µs,
near-field pair force 22.139 µs, one mobility pass 151.673 µs, one contact pass
85.646 µs, whole physical step 1514.281 µs. With two mobility and twelve contact passes,
contact remains the largest recurring expense. These isolated timings are indicative,
not asserted additive wall-time accounting. Sample-download timing includes waiting
for preceding GPU work; it is not pure PCIe transfer cost. Cold compilation/setup
are separately recorded. Small core logs show about 1.25–1.27 s cold F2 compilation
for the exercised variants; cached compilation times in later summaries are much smaller.

The regression stays below one second and within the same order as F1, but the
2.06× total cost is reported rather than hidden. The historical 6.3× batching speedup
is **not reclassified as an F2 measurement**. Batch independence was checked on tiny
samples; the architecture/config are preserved. No new long run or batch campaign
was required.

## 8. Deliverables, scope and next evidence needed

Four English PNG/PDF groups with source arrays are in `results/fast_pzt/f2/figures`:
`Fig_F2_1_near_field`, `Fig_F2_2_contact_sensitivity`,
`Fig_F2_3_PZT_morphodynamics`, `Fig_F2_4_3D_cluster_shapes`.
The MP4 is `visualization/F2_PZT_mature_cluster_actual_radius_60fps.mp4`:
1112 decoded frames, 60 fps, 1200×900, 0–222.2 ms actual time, 18.533 s playback.
Each frame is an actual stored state. The common precursor and its true exact-state
continuation are indexed in `STATE_PROVENANCE.json`; no interpolated state, changed
sphere radius, cluster magnification, ID subsampling or moving camera is used.

The response domain remains radius 1.779–1.781 µm, frequency 2.499–2.501 MHz and
pressure 0.1–40000 Pa; this campaign uses 1.78 µm, 2.5 MHz only. The basis supports
2.5 MHz only. Input SHA-256 is in `INPUT_MANIFEST.json`; both inputs are unchanged.
No flow, walls, polydispersity, shape deformation, added mass/history, collective
radial ODE, clot response or calibrated hardware acoustics is claimed.

Next evidence should target near-contact force/phase validity, angularly consistent
dense hydrodynamics/contact, and a reproducible morphology envelope before interpreting
micro-packing or transient axes quantitatively. A matched theory/experiment is needed
for calibration. This delivery does not start those further tasks.

Reproduction commands and data semantics are in
[`scripts/fast_pzt_f2/README.md`](../scripts/fast_pzt_f2/README.md).
`ARTIFACT_MANIFEST.json` indexes closed files by SHA-256 and local/server location.
Necessary small raw trajectories, IDs, logs, figures and movie are archived in Git;
unchanged large input files remain indexed externally. Failed contact/time-step runs
are retained and excluded explicitly from selected-backend quality claims.
