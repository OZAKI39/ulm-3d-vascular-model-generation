# F2 implementation and frozen experiment plan

Base archive `6de76d508ea088e137069abf058b86b7f6d47a1c`; historical F1 scientific
execution remains `9829156893339efff713483f3125848b84adbc98` (matched controls
`3721bebd597164262691af766f723d245e043116`). New branch is
`dev/fast-pzt-aggregation-f2-20261008`. All six F1 simulator modules remain byte
identical. F2 subclasses the engine and overrides only the pair/contact stages;
new results and supervision use a separate F2 root. No active F1 job is resumed.

## Candidate hypotheses, fixed before trajectories

Near-field candidate: a **pairwise fundamental monopole-feedback closure**, not
Doinikov's multipole translational equations or Fan's reduced cluster radial ODE.
For positive-exponential phasors, p_sc,j = -K V_j, K = rho*omega^2/(4*pi*d).
Using the table's local fundamental secant compliance chi_i = V_i0/P_i,
solve the two scalar complex equations V_i + K chi_i V_j = V_i0 and
V_j + K chi_j V_i = V_j0 analytically. Use the resulting signed fundamental
correlation only for the pair force; retain harmonics 2–8 unchanged. Blend the
correlation with the uncoupled value by a cubic smoothstep: full feedback at
contact, zero at d = 4(Ri+Rj). Retain F1's geometric regularizer in both models.
There is no fitted attraction multiplier or sign clipping. Algebraic response
phase changes can change a general pair's signed correlation. In-phase and
counterphase identical pairs must retain attraction/repulsion, respectively.
At zero external pressure, use the table's lowest valid-pressure compliance as
a declared linear-limit approximation; all-OFF still gives exactly zero response.
Reject singular closure denominators and coupled effective pressure outside the
original response domain. This is a force-only pair closure: it does not create a
consistent collective bubble response or modify Primary force. Its deficiencies
are part of the sensitivity experiment; physics remains EXPLORATORY.

Contact candidate: freeze local dimensionless mobility tensors
M_i = [I + sum_j c_ij n_ij n_ij^T]^-1 at the predicted state, using exactly F1's
normal resistance coefficient. For each overlapping pair use
H_ij = (M_i+M_j)/2 and equal/opposite corrections
δ_i = (h_overlap/2) H_ij n/(n^T H_ij n), δ_j = -δ_i, with F1's 0.5 relaxation
and 12 Jacobi passes. This solves the pair's linearized constraint in a symmetric
relative-mobility metric. It preserves center of mass, responds to local crowding,
and has positive separating normal work. Anisotropic response is not friction.
The sum of pair corrections cannot reduce global squared Rg in a contact pass.
It is not the full coupled Stokes contact metric. No stiffness fitting, adhesion,
cluster locking, velocity cap or global QP is added.

## Small, ordered comparisons

A: distance/force curve, analytic 2x2 closure versus Warp, same/counterphase signs,
then two actual two-bubble trajectories. B: one fixed-seed 32-bubble scene,
F1/F1 versus F2-force/F1-contact, then F2-force/F2-contact; one matched OFF case.
F1's archived trajectory supplies the baseline, with a tiny F1/F2-delegation
check instead of repeating its campaign. If a candidate fails 20 us quality,
try only 10 us then, if necessary, 5 us on this representative scene; matched
comparisons use identical dt. If still invalid, preserve it and retain the
qualified F1 component. No force parameters are retuned to match morphology.

C: form a real-basis cluster within 160 ms. Maturity requires >=90% membership,
unchanged largest-member set and connected membership for 10 ms, Rg variation
<=2% across that window, and accepted overlap <=0.01R. Validate topology against
all accepted-step edge changes, not just movie frames. Save the first qualifying
actual state. From it run unchanged, phase-only and amplitude/enable controls
for 120 ms each. Select controls by a short field-only diagnostic before these
trajectories: comparable cloud RMS pressure (+/-10%), maximize normalized
complex-gradient contrast. Phase-only keeps source amplitudes/enables exactly
fixed. At most three phase and three amplitude/enable candidates are diagnosed.
No second control search after observing morphology.

Use one unchanged-control repeat for floating-point sensitivity and, where
qualified, one unchanged-control branch per alternate force/contact component
from the same mature state, plus one two-versus-four local mobility iteration check. These bounded, single-factor checks compare control
effects with numerical/model sensitivity without a factorial campaign.

D: one 200-bubble/5-ms regression with the selected F2 backend; reuse exact F1
seed/cloud/pressure. Decompose pair/mobility/contact cost if substantially slower.
Preserve batching and check only a tiny independent two-sample state.

## Outputs and interpretation

20 us, one actual substep, two mobility iterations and 12 contact passes remain
nominal. Preserve F1 overlap/correction/finite/domain thresholds. Freeze new
observation thresholds in `configs/fast_pzt/f2/protocol.json` before execution.
Keep raw connectivity events; classify persistence separately (2 ms dwell,
right-censor unobserved end intervals). Store initial fixed-member shape and
dynamic component shape, rigid Kabsch rotation, aligned internal RMSD, Rg,
eigenvalues, principal-axis reliability and adjacency turnover. Centroid change
alone is never shape evidence. All plots use actual positions/identities/radii;
the movie joins only an actual precursor and its exact-state continuation, with
an explicit provenance map and no interpolated or fictitious states.

Four English figures, one small 3D MP4, model/contact/control/performance JSONs,
separate gates, report, source/input hashes and reproducible commands are delivered.
The F1 result archive is read-only. There is no claim of quantitative calibration.

References: https://doi.org/10.1103/PhysRevE.92.043001 (publisher abstract inspected;
full PDF unavailable through tool), https://doi.org/10.1016/j.ultsonch.2019.104790
(scope reference, not implemented), https://mmacklin.com/xpbd.pdf (constraint
metric concept; inertial mass weights are not imported into the overdamped model).

## Contact candidate disposition after the first bounded comparison

The tensor pair metric fails the original overlap budget at 20/10 us and passes
at 5 us. A subsequent algebraic check exposes a nonzero angular moment of its
position correction (equal/opposite anisotropic shifts need not be central).
This invariant defect is independent of the desired morphology. Preserve the
candidate and its data for sensitivity analysis, but reject it for causal PZT
experiments. Adopt the new near-field force with F1 central PBD at the first
qualified 10 us step. Compare near fields at 10 us and contacts at 5 us, each with
matched dt and unchanged quality criteria. Do not implement another contact
variant or tune the force to recover 20 us. One contact-diagnostic mature branch
may be retained as explicitly unqualified evidence, not in the control effect
claim or qualified dataset. The seven-branch limit is unchanged.

Field-only eligibility repair: the three initial phase patterns changed ROI RMS
pressure by 16.7–27.5%, exceeding the unchanged 10% matching bound. Preserve those
records. Before the phase-only trajectory, set Y+ to +30 degrees and solve the
single trigonometric RMS equation for Y−, choosing the root nearest zero. Source
amplitudes/enables remain exactly unchanged. This adds one analytic field control,
not a morphology-driven search, and keeps the same mature state and observation
criteria. Record it in a separate run rather than rewriting the original selection.
