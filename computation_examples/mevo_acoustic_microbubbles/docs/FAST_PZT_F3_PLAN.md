# Frozen F3 scope and ordered experiment plan

Base archive 6655867b772484a82ea9a7c2afd491da02bdd378; historical F2 execution commits
remain attached to their datasets. Branch dev/fast-pzt-aggregation-f3-20261008.
Physics, geometry, input table/basis and tolerances are frozen. No old job resumes.

A first: evaluate the exact F2 production force/mobility prefix at the saved 102.2 ms
state, without contact or position advance. Separate Primary, Secondary, total,
control-induced centered differences, free Stokes and local/balanced mobility velocity.
Fit affine gradients with rank/condition checks. F_scale is Stokes drag for 1 um/s,
not a fitted force floor. Diagnose occupied grid cells, actual complex fields,
spatial gradients, local finite differences of interpolated pressure and consequent
Primary-force differences. The stored-gradient interpolation is retained even if its
local consistency is imperfect; report this uncertainty, do not silently repair physics.

Screen six configurations only: unchanged, the F2 phase-only control, one analytically
RMS-matched Z-pair phase control (Z+ = 30 degrees; solve Z−, smallest absolute root),
and the three previously screened F2 amplitude/enable vectors. Match/report both actual
bubble-position and padded cluster-ROI RMS pressure. No hardware amplitude limit is
stored in the fast interface; use the conservative maximum previously field-screened
F2 source amplitude, 291.237687721856 Pa, as the F3 campaign bound. No hardware rating
is inferred. Select one phase and one amplitude candidate by predeclared instantaneous
mobility deviatoric strain contrast, after domain/amplitude/pressure eligibility. Keep
raw forces and both absolute/relative values; the diagnostic ratios are not calibrated
shape thresholds. The ROI is the actual bubble bounding box padded by 2R, 5³ points.

B next: preserve fast_atomic exactly. New stable_reduction stores every unique F2
pair contribution once and sums partners in ascending ID order; component velocity
balance sums members in ascending ID order. Contact/mobility iterations and the
integer component union graph remain unchanged. No compensated/new physics arithmetic
is added. Diagnose eight fixed-state repetitions per mode and reduction ablations,
then three 120 ms unchanged-control repeats per primary mode. Check fixed-state
NumPy/GPU force equivalence, deterministic component labels, graph versus direct
small-step replay, exact control clock and numerical quality. Attribute only what
these experiments resolve; do not assume every difference is caused by atomics.

C/D: after screening and reproducibility, run only the selected phase and amplitude
continuations from exactly the same archived state for 120 ms in the qualified mode;
reuse its unchanged repeat as reference. Two predeclared zero-mean initial perturbations
(1 pm and 1 nm RMS, same fixed direction seed) separately check true initial-condition
sensitivity if quality permits; do not project an invalid starting state. All 32 IDs
remain free. Use F2's absolute shape thresholds plus 3× the maximum pairwise repeat
difference. A sustained signal requires 10 ms; topology dwell remains 2 ms. Endpoint,
short-window, all-time and transient results are separate. Preserve the inherited F2
near-field/mobility sensitivity envelope as a conservative, historical model comparison;
do not rerun or change those physics. Numerical repeatability is not physical accuracy.

One N200/5ms regression per execution mode quantifies their separate costs. No larger
campaign. Four English figures use actual spheres and fixed camera/scale, with explicit
force-vector scaling. A short actual-frame comparison movie may illustrate a negative
or non-identifiable result. Closed raw trajectories, execution identities, failed runs,
logs, hashes, report and module gates are archived; old F1/F2 files are read-only.

The Z-pair analytic RMS equation has no feasible root on the actual small ROI in
the initial field attempt. That candidate slot is marked ineligible, with Z−=0
retained only for field diagnostics. No extra control replaces it; the six-slot
budget and original F2 phase candidate are retained. Failed startup/field attempts
remain logged. No dynamic result was observed before this disposition.
